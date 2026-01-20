use crate::{Real, Complex};
use super::{Constraint, Variable, Assignment, ConstraintLearning, Clause};
use batsat::{Solver, SolverOpts, BasicCallbacks, lbool, SolverInterface};

#[derive(Debug, Clone)]
pub struct MaybeData<C: Complex> {
	pub area: C::Real,
	pub max_mod2: C::Real,
}

impl<C: Complex> MaybeData<C> {
	fn merge(from: Vec<Self>) -> Self {
		let mut total_area = C::Real::zero();
		let mut total_max_mod2 = C::Real::zero();
		for data in from {
			total_area = total_area + data.area;
			total_max_mod2 = total_max_mod2 + data.max_mod2;
		}
		MaybeData {
			area: total_area,
			max_mod2: total_max_mod2,
		}
	}
}

#[derive(Debug, Clone)]
pub enum QsatResult<C: Complex> {
	Sat,
	Unsat,
	Maybe(MaybeData<C>),
}

#[derive(Debug, Clone)]
pub(super) enum TheoryResult<C: Complex> {
	Maybe(Vec<MaybeData<C>>),
	Conflict(Vec<Clause>),
}

#[derive(Debug)]
pub struct QsatInstance<C: Complex> {
	num_qubits: u16,
	max_depth: u16,
	promise_min_sat: Option<C::Real>,
	constraints: Vec<Constraint<C>>,
}

impl <C: Complex> QsatInstance<C> {
	pub fn new(num_qubits: u16, max_depth: u16,) -> Self {
		QsatInstance {
			num_qubits,
			max_depth,
			constraints: Vec::new(),
			promise_min_sat: None,
		}
	}

	pub fn add_constraint(&mut self, support: Vec<u16>, projector: Vec<C>) {
		let max_qubit = *support.iter().max().unwrap();
		assert!(max_qubit < self.num_qubits, "qubit index out of range");
		assert!(projector.len() == 2usize.pow(support.len() as u32),
			"projector size does not match the number of qubits");
		self.constraints.push(Constraint::new(support, projector));
	}

	pub fn set_promise(&mut self, promise: C::Real) {
		self.promise_min_sat = Some(promise);
	}

	pub fn check_sat(&self, result: QsatResult<C>) -> QsatResult<C> {
		match result {
			QsatResult::Maybe(data) => {
				if let Some(promise) = &self.promise_min_sat
				&& C::Real::int(self.num_qubits.into()) * data.max_mod2.clone() < *promise {
					QsatResult::Sat
				} else {
					QsatResult::Maybe(data)
				}
			}
			_ => { result }
		}
	}

	#[allow(dead_code)]
	fn constraints_involving(&self, qubit: u16) -> Vec<usize> {
		let mut result = Vec::new();
		for (pos, constraint) in self.constraints.iter().enumerate() {
			if constraint.support().contains(&qubit) {
				result.push(pos);
			}
		}
		result
	}

	fn num_vars(&self) -> u32 {
		(self.num_qubits as u32) * (2 * (self.max_depth as u32) - 1)
	}

	fn assignment(&self) -> Assignment {
		Assignment::new(self.num_qubits, self.max_depth)
	}

	fn learn(&self, assignment: Assignment) -> TheoryResult<C> {
		let mut maybes = Vec::new();
		let mut conflicts = Vec::new();
		for c in &self.constraints {
			let result = assignment.find_conflict(c, true);
			match result {
				ConstraintLearning::Maybe(data) => {maybes.push(data);},
				ConstraintLearning::Conflict(clause) => {conflicts.push(clause);},
			}
		}
				
		if !conflicts.is_empty() {
			TheoryResult::Conflict(conflicts)
		} else {
			TheoryResult::Maybe(maybes)
		}
	}

}

pub fn solve<C: Complex>(instance: QsatInstance<C>) -> QsatResult<C> {
	let mut solver = Solver::new(SolverOpts::default(), BasicCallbacks::new());
	let mut assignment = instance.assignment();
	for i in 0..instance.num_vars() {
		let v = solver.new_var(lbool::TRUE, true);
		assert!(v.idx() == i);
	}
	debug_assert!(solver.num_vars() == instance.num_vars());


	let result = loop {
		let mut temp_assignment = assignment.clone();
		temp_assignment.fill();
		match instance.learn(temp_assignment) {
			TheoryResult::Maybe(maybes) => {
				let data = MaybeData::merge(maybes);
				break QsatResult::Maybe(data);
			}
			TheoryResult::Conflict(clauses) => {
				for clause in clauses {
					let mut sat_clause = clause.to_sat_clause(instance.num_qubits);
					solver.add_clause_reuse(&mut sat_clause);
				}
			}
		}

		if solver.solve_limited(&[]) == lbool::FALSE {
			break QsatResult::Unsat;
		}
		assignment.update_from_solver(&solver);
	};

	instance.check_sat(result)
}