use crate::checker::{check_constraint, ConstraintResult, AngleFraction};
use crate::{Complex, MaybeData};

#[derive(Copy, Clone,Debug, PartialEq, Eq)]
pub enum AngleId{
	Theta(u16),
	Phi(u16),
}

impl AngleId {
	pub fn set<I>(self, vec: &mut Vec<(I, I)>) -> &mut I {
		match self {
			AngleId::Theta(q) => &mut vec[q as usize].0,
			AngleId::Phi(q) => &mut vec[q as usize].1,
		}
	}
	pub fn get <I: Copy>(self, vec: &Vec<(I, I)>) -> I {
		match self {
			AngleId::Theta(q) => vec[q as usize].0,
			AngleId::Phi(q) => vec[q as usize].1,
		}
	}
}


#[derive(Copy, Clone, Debug)]
pub(super) struct Variable {
	qubit: u16,
	is_phi: bool,
	depth: u16,
}

impl Variable {
	pub(super) fn index(&self, num_qubits: u16) -> u32 {
		let depth = self.depth as u32;
		let nqubits = num_qubits as u32;
		let q = self.qubit as u32;
		let is_theta = !self.is_phi as u32;
		let result = depth * (nqubits * 2) + nqubits * is_theta + q;
		result
	}
	pub(super) fn from_index(index: usize, num_qubits: u16) -> Self {
		let depth = (index as u32 / (num_qubits as u32 * 2)) as u16;
		let rem = index as u32 % (num_qubits as u32 * 2);
		if rem >= num_qubits as u32 {
			let qubit = (rem - num_qubits as u32) as u16;
			return Variable { qubit, is_phi: false, depth };
		} else {
			let qubit = rem as u16;
			return Variable { qubit, is_phi: true, depth };
		}
	}
}

#[derive(Debug, Clone)]
pub(super) struct Clause(pub Vec<(Variable, bool)>);

impl Clause {
	pub fn empty() -> Self {
		Clause(Vec::new())
	}
	pub fn push(&mut self, var: Variable, value: bool) {
		self.0.push((var, value));
	}
}

#[derive(Debug)]
pub(super) enum ConstraintLearning<C: Complex> {
	Maybe(MaybeData<C>),
	Conflict(Clause),
}

impl<C: Complex> ConstraintLearning<C> {
	#[cfg(test)]
	pub fn is_conflict(&self) -> bool {
		matches!(self, ConstraintLearning::Conflict(_))
	}
}

#[derive(Clone, Debug)]
pub(super) struct Assignment {
	phis: Vec<Vec<Option<bool>>>,
	thetas: Vec<Vec<Option<bool>>>,
} // angle[nqubit][depth]

impl Assignment {
	pub fn new(nqubits: u16, depth: u16) -> Self {
		let nqubits = nqubits as usize;
		let depth = depth as usize;
		Assignment {
			phis: vec![vec![None; depth]; nqubits],
			thetas: vec![vec![None; depth - 1]; nqubits],
		}
	}

	pub(super) fn num_qubits(&self) -> u16 {
		self.phis.len() as u16
	}

	fn get_qubit_angles(&self, qubit: u16, depth: Option<(u16, u16)>) -> (AngleFraction, AngleFraction) {
		let theta = get_angle_fraction(
			&self.thetas[qubit as usize], true, depth.map(|d| d.0));
		let phi = get_angle_fraction(
			&self.phis[qubit as usize], false, depth.map(|d| d.1));
		(theta, phi)
	}

	fn get_angles(&self, qubits: &Vec<u16>) -> Vec<(AngleFraction, AngleFraction)> {
		qubits.iter().map(|q| self.get_qubit_angles(*q, None)).collect()
	}


	fn get_angles_depth(&self, qubits: &Vec<u16>, depths: &Vec<(u16, u16)>)
	-> Vec<(AngleFraction, AngleFraction)> {
		qubits.iter().zip(depths.iter())
			.map(|(q, d)| self.get_qubit_angles(*q, Some(*d)))
			.collect()
	}

	pub fn fill(&mut self) -> bool {
		let mut changed = false;
		for item in [&mut self.phis, &mut self.thetas].iter_mut().flat_map(|v| v.iter_mut()).flatten(){
			if item.is_none() {
				*item = Some(false);
				changed = true;
			}
		}
		changed
	}

	pub fn assign(&mut self, var: Variable, value: Option<bool>) {
		let vec = if var.is_phi {
			&mut self.phis[var.qubit as usize]
		} else {
			&mut self.thetas[var.qubit as usize]
		};

		if vec.len() <= var.depth as usize {
			vec.resize(var.depth as usize + 1, None);
		}
		vec[var.depth as usize] = value;
	}

	pub fn get(&self, var: Variable) -> Option<bool> {
		let vec = if var.is_phi {
			&self.phis[var.qubit as usize]
		} else {
			&self.thetas[var.qubit as usize]
		};
		if vec.len() <= var.depth as usize {
			None
		} else {
			vec[var.depth as usize]
		}
	}

	fn depths(&self, qubit: u16) -> (u16, u16) {
		let phi_vec = &self.phis[qubit as usize];
		let theta_vec = &self.thetas[qubit as usize];
		let phi_depth = phi_vec.iter()
			.position(|x| x.is_none()).unwrap_or(phi_vec.len());
		let theta_depth = theta_vec.iter()
			.position(|x| x.is_none()).unwrap_or(theta_vec.len());
		(theta_depth as u16, phi_depth as u16)
	}

	pub fn find_conflict<C: Complex>(&self, constraint: &Constraint<C>, thorough: bool)
	-> ConstraintLearning<C> {
		if let ConstraintResult::Maybe(data) = constraint.check(self.get_angles(&constraint.support)) {
			return ConstraintLearning::Maybe(data);
		}
		let depths = if thorough {
			self.find_conflict_depths::<C>(&constraint)
		} else {
			constraint.support.iter().map(|q| self.depths(*q)).collect()
		};

		let mut clause = Clause::empty();
		let support  = &constraint.support;
		for (q, depths) in support.iter().zip(depths.iter()) {
			let q = *q;
			let (theta_depth, phi_depth) = depths;
			if *theta_depth > 0 {
				let var = Variable { qubit: q, is_phi: false, depth: *theta_depth - 1 };
				clause.push(var, !self.get(var).unwrap())
			}
			if *phi_depth > 0 {
				let var = Variable { qubit: q, is_phi: true, depth: *phi_depth - 1 };
				clause.push(var, !self.get(var).unwrap())
			}
		}
		ConstraintLearning::Conflict(clause)
	}

	fn find_conflict_depths<C: Complex>(&self, constraint: &Constraint<C>) -> Vec<(u16, u16)> {

		let support = &constraint.support;
		let mut depths: Vec<(u16, u16)> = Vec::with_capacity(support.len());
		let mut can_move = Vec::with_capacity(support.len() * 2);

		for q in support {
			can_move.push(AngleId::Theta(*q));
			can_move.push(AngleId::Phi(*q));
			let depth = self.depths(*q);
			depths.push(depth);
		}
		can_move.retain(|item| item.get(&depths) > 0);

		while !can_move.is_empty() {
			let mut new_depths = depths.clone();
			for item in &can_move {
				*item.set(&mut new_depths) -=  1;
			}
			if constraint.check(
				self.get_angles_depth(&constraint.support, &new_depths)
			).is_unsat() {
				depths = new_depths;
				can_move.retain(|x| x.get(&depths) > 0);
				continue;
			} else {
				can_move.sort_unstable_by(|a, b| b.get(&depths).cmp(&a.get(&depths)));
				let mut new_can_move = Vec::new();
				for item in can_move {
					let mut new_depths = depths.clone();
					*item.set(&mut new_depths) -= 1;
					if constraint.check(
						self.get_angles_depth(&constraint.support, &new_depths)
					).is_unsat() {
						if *item.set(&mut new_depths) > 0 {
							new_can_move.push(item);
						}
						depths = new_depths;
					}
				}
				can_move = new_can_move;
			}
		}

		depths
	}
}

fn get_angle_fraction(path: &Vec<Option<bool>>, first_left: bool, max_depth: Option<u16>)
-> AngleFraction {
	let mut numerator = 0usize;
	let mut denominator = if first_left {2} else { 1 } as usize;
	for (depth, bit) in path.iter().enumerate() {
		if let Some(max_depth) = max_depth {
			if depth as u16 >= max_depth {
				break;
			}
		}
		if let Some(bit) = bit {
			numerator *= 2;
			denominator *= 2;
			if *bit {
				numerator += 1;
			}
		} else { break; }
	}
	AngleFraction::new(numerator, denominator)
}

#[derive(Debug, Clone)]
pub(super) struct Constraint<C: Complex> {
	support: Vec<u16>,
	projector: Vec<C>,
}

impl<C: Complex> Constraint<C> {
	fn check(&self, angles: Vec<(AngleFraction, AngleFraction)>) -> ConstraintResult<C> {
		let result = check_constraint::<C>(&self.projector, angles.clone());
		result
	}
	pub(super) fn new(support: Vec<u16>, projector: Vec<C>) -> Self {
		Constraint { support, projector }
	}
	pub(super) fn support(&self) -> &Vec<u16> {
		&self.support
	}
}

#[cfg(test)]
mod test {
	use super::*;
	use crate::projectors::computational;
	use crate::numbers::ComplexFloat;
	#[test]
	fn test_clause_learning() {
		let projector = computational::<ComplexFloat>(0, 3);
		let constraint = Constraint {
			support: vec![0,1,2],
			projector,
		};

		let mut assignment = Assignment::new(3, 7);
		assignment.fill();

		let learning = assignment.find_conflict::<ComplexFloat>(&constraint, true);
		assert!(learning.is_conflict());
		let clause = match learning {
			ConstraintLearning::Conflict(c) => c,
			_ => panic!("expected conflict"),
		};

		let mut assignment = Assignment::new(3, 7);
		for (var, value) in clause.0 {
			assignment.assign(var, Some(!value))
		}

		assert!(assignment.find_conflict(&constraint, true).is_conflict());
	}
}