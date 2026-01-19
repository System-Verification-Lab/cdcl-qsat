use oxiz_sat::{Solver, SolverResult, Var, Lit, TheoryCallback, TheoryCheckResult};
use smallvec::SmallVec;
use crate::checker::{check_constraint, ConstraintResult, Fraction};
use crate::{Complex, Real};

#[derive(Copy, Clone, Debug)]
struct Variable {
	qubit: u16,
	is_phi: bool,
	depth: u8,
}

impl Variable {
	fn index(&self, num_qubits: u16) -> u32 {
		self.depth as u32 * (num_qubits as u32 * 2) + (self.qubit as u32 * 2) + (self.is_phi as u32)
	}

	fn from_index(index: usize, num_qubits: u16) -> Self {
		let depth = (index as u32 / (num_qubits as u32 * 2)) as u8;
		let rem = index as u32 % (num_qubits as u32 * 2);
		let qubit = (rem / 2) as u16;
		let is_phi = (rem % 2) == 1;
		Variable { qubit, is_phi, depth }
	}
	fn qubit(&self) -> u16 {
		self.qubit
	}
	fn new(qubit: u16, is_phi: bool, depth: u8) -> Self {
		Variable { qubit, is_phi, depth }
	}
}

#[derive(Debug, Clone)]
pub enum QsatResult<C: Complex> {
	Sat,
	Unsat,
	Maybe {sum_area: C::Real, sum_max_mod2: C::Real},
}

#[derive(Clone, Debug)]
struct QubitAssignment {
	phi: Vec<Option<bool>>,
	theta: Vec<Option<bool>>,
}

impl QubitAssignment {
	pub fn new() -> Self {
		QubitAssignment {
			phi: Vec::new(),
			theta: Vec::new(),
		}
	}

	pub fn theta(&self) -> Fraction {
		get_interval_start(&self.theta, false)
	}

	pub fn phi(&self) -> Fraction {
		get_interval_start(&self.phi, true)
	}

	pub fn assign(&mut self, is_phi: bool, depth: u16, value: Option<bool>) {
		let vec = if is_phi { &mut self.phi } else { &mut self.theta };
		if vec.len() <= depth as usize {
			vec.resize(depth as usize + 1, None);
		}
		vec[depth as usize] = value;
	}

	pub fn negate(&self, qubit: u16, num_qubits: u16) -> SmallVec<[Lit; 8]> {
		let mut clause = SmallVec::<[Lit; 8]>::new();
		build_clause(&mut clause, &self.theta, qubit, num_qubits, false);
		build_clause(&mut clause, &self.phi, qubit, num_qubits, true);
		clause
	}
	pub fn get_unassigned_var(&self, qubit: u16, num_qubits: u16) -> Option<Var> {
		for (depth, bit) in self.theta.iter().enumerate() {
			if bit.is_none() {
				return Some(Var::new(
					Variable::new(qubit, false, depth as u8).index(num_qubits)
				));
			}
		}
		for (depth, bit) in self.phi.iter().enumerate() {
			if bit.is_none() {
				return Some(Var::new(
					Variable::new(qubit, true, depth as u8).index(num_qubits)
				));
			}
		}
		None
	}
}

fn get_interval_start(path: &Vec<Option<bool>>, first_left: bool) -> Fraction {
	let mut numerator = 0usize;
	let mut denominator = if first_left {2} else { 1 } as usize;
	for bit in path {
		if let Some(bit) = bit {
			numerator *= 2;
			denominator *= 2;
			if *bit {
				numerator += 1;
			}
		} else { break; }
	}
	Fraction::new(numerator, denominator)
}

fn build_clause(
	clause: &mut SmallVec<[Lit; 8]>, bits: &Vec<Option<bool>>,
	qubit: u16, num_qubits: u16, is_phi: bool
) {
	for (depth, bit) in bits.iter().enumerate() {
		if let Some(bit) = bit {
			let var = Var::new(
				Variable::new(qubit, is_phi, depth as u8).index(num_qubits)
			);
			let lit = if *bit { Lit::neg(var) } else { Lit::pos(var) };
			clause.push(lit);
		} else { break; }
	}
}

#[derive(Debug)]
struct Constraint<C: Complex> {
	qubits: Vec<u16>,
	projector: Vec<C>,
}

#[derive(Debug)]
pub struct QsatSolver<C: Complex> {
	num_qubits: u16,
	max_depth: u16,
	assignments: Vec<QubitAssignment>,
	state: Vec<Option<ConstraintResult<C>>>,
	assigned_on_levels: Vec<Vec<Variable>>,
	promise_min_sat: Option<C::Real>,
	constraints: Vec<Constraint<C>>,
	level: u32,
}

impl <C: Complex> QsatSolver<C> {
	pub fn new(num_qubits: u16, max_depth: u16,) -> Self {
		let total_vars = (num_qubits as u64) * 2 * (max_depth as u64);
		assert!(total_vars <= u32::MAX as u64, "too many variables for u32 indices");
		QsatSolver {
			num_qubits,
			max_depth,
			assignments: vec![QubitAssignment::new(); num_qubits as usize],
			assigned_on_levels: Vec::new(),
			constraints: Vec::new(),
			state: Vec::new(),
			level: 0,
			promise_min_sat: None,
		}
	}
	pub fn add_constraint(mut self, qubits: Vec<u16>, projector: Vec<C>) -> Self {
		let max_qubit = *qubits.iter().max().unwrap();
		assert!(max_qubit < self.num_qubits, "qubit index out of range");
		assert!(projector.len() == 2usize.pow(qubits.len() as u32),
			"projector size does not match the number of qubits");
		self.constraints.push(Constraint { qubits, projector });
		self.state.push(None);
		self
	}

	pub fn solve(mut self) -> QsatResult<C> {
		let mut solver = Solver::new();
		let nvars = (self.num_qubits as usize) * 2 * (self.max_depth as usize);
		solver.ensure_vars(nvars);
		for i in 0..nvars {
			let var = Var::new(i as u32);
			let lit_pos = Lit::pos(var);
			let lit_neg = Lit::neg(var);
			solver.add_clause([lit_pos, lit_neg]);
		}
		let result = solver.solve_with_theory(&mut self);
		dbg!(solver.trail());
		dbg!(&self.assignments);
		match result {
			SolverResult::Unsat => QsatResult::Unsat,
			_ => {
				let mut sum_max_mod2 = C::Real::zero();
				let mut sum_area = C::Real::zero();
				for constraint_result in &self.state {
					if let Some(ConstraintResult::Maybe { area, max_mod2, .. }) = constraint_result {
						sum_area = sum_area + area.clone();
						sum_max_mod2 = sum_max_mod2 + max_mod2.clone();
					}
				}
				if let Some(promise) = &self.promise_min_sat
				&& C::Real::int(self.state.len()) * sum_max_mod2.clone() < *promise {
					QsatResult::Sat
				} else {
					QsatResult::Maybe {sum_area, sum_max_mod2}
				}	
			}
		}
	}

	fn check_constraint(&mut self, pos: usize) -> &ConstraintResult<C> {
		let constraint = &self.constraints[pos];
		let mut angle_limits = Vec::new();
		for qubit in &constraint.qubits {
			let assignment = &self.assignments[*qubit as usize];
			let theta = assignment.theta();
			let phi = assignment.phi();
			angle_limits.push((theta, phi));
		}
		let result = check_constraint::<C>(&constraint.projector, angle_limits);
		self.state[pos] = Some(result);
		&self.state[pos].as_ref().unwrap()
	}

	fn unassign(&mut self, var: Variable) {
		let assignment = &mut self.assignments[var.qubit as usize];
		assignment.assign(var.is_phi, var.depth as u16, None);
		for pos in self.constraints_involving(var.qubit) {
			self.state[pos] = None;
		}
	}

	fn assign(&mut self, var: Variable, value: bool) {
		let assignment = &mut self.assignments[var.qubit as usize];
		assignment.assign(var.is_phi, var.depth as u16, Some(value));
		self.assigned_on_levels.last_mut().unwrap().push(var.clone());
	}

	fn constraints_involving(&self, qubit: u16) -> Vec<usize> {
		let mut result = Vec::new();
		for (pos, constraint) in self.constraints.iter().enumerate() {
			if constraint.qubits.contains(&qubit) {
				result.push(pos);
			}
		}
		result
	}

	fn check_constraints(&mut self, which: impl Iterator<Item = usize>, recheck: bool)
	-> TheoryCheckResult {
		let mut clause: Option<SmallVec<[Lit; 8]>> = None;
		for pos in which {
			let result = if recheck || self.state[pos].is_none() {
				self.check_constraint(pos)
			} else {
				self.state[pos].as_ref().unwrap()
			};
			if result.is_unsat() && clause.is_none() {
				let mut new_clause = SmallVec::<[Lit; 8]>::new();
				for qubit in self.constraints[pos].qubits.iter() {
					new_clause.extend(
						self.assignments[*qubit as usize].negate(*qubit, self.num_qubits)
					);
				}
				clause = Some(new_clause);
			}
		}


		if let Some(clause) = clause {
			TheoryCheckResult::Conflict(clause)
		} else if !recheck {
			let mut clause = None;
			for (qubit, assignment) in self.assignments.iter().enumerate() {
				if let Some(unassigned) = assignment.get_unassigned_var(qubit as u16, self.num_qubits) {
					println!("Found unassigned variable: {:?}", unassigned);
					let mut new_clause = SmallVec::<[Lit; 8]>::new();
					new_clause.push(Lit::pos(unassigned));
					//new_clause.push(Lit::neg(unassigned));
					clause = Some(new_clause);
					break;
				}
			}
			if let Some(clause) = clause {
				println!("Adding clause to force assignment: {:?}", clause);
				TheoryCheckResult::Conflict(clause)
			} else { TheoryCheckResult::Sat}
		} else {
			TheoryCheckResult::Sat
		}
	}
}

impl<C: Complex> TheoryCallback for QsatSolver<C> {
	fn on_assignment(&mut self, lit: Lit) -> TheoryCheckResult {
		let index = lit.var().index();
		println!("Assigning variable index {}", index);
		let value = lit.is_pos();
		let var = Variable::from_index(index, self.num_qubits);
		let qubit = var.qubit();
		self.assign(var, value);
		let constraints = self.constraints_involving(qubit);
		self.check_constraints(constraints.into_iter(), true)
	}	

	fn final_check(&mut self) -> TheoryCheckResult {
		println!("Final check!");
		self.check_constraints((0..self.constraints.len()).into_iter(), false)
	}

	fn on_backtrack(&mut self, level: u32) {
		print!("Backtracking to level {}\n", level);
		for _ in (level..self.level - 1).rev() {
			let vars_to_unassign = self.assigned_on_levels.pop().unwrap();
			for var in vars_to_unassign {
				self.unassign(var);
			}
		}
		self.level = level;
	}

	fn on_new_level(&mut self, level: u32) {
		print!("Entering level {}\n", level);
		self.assigned_on_levels.push(Vec::new());
		self.level = level;
	}

}
