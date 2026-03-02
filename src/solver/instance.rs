use crate::{Complex, Real, solver::instance};
use super::{Constraint, Variable, Assignment, ConstraintLearning, Clause};
use cpu_time::{ProcessTime};
use rustsat::solvers::{Solve, SolveStats, SolverResult};
use rustsat::types::{Lit, TernaryVal, Var};
use rustsat_batsat::BasicSolver as BatSat;
use rustsat_cadical::CaDiCaL;
use rustsat_glucose::core::Glucose;
use std::time::Duration;
use std::marker::PhantomData;


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
pub struct QsatDiagnosis<C: Complex> {
	pub result: QsatResult<C>,
	pub execution_time: Duration,
	pub num_vars: usize,
	pub num_calls_solver: usize,
	pub num_clauses: usize,
	pub num_cdcl_constraints: Option<usize>,
}

pub trait AbstractSatSolver: Solve + Default + SolveStats {
	fn clause_count(&self) -> usize;
	fn label(&self) -> &'static str;
	fn cdcl_constraints(&self) -> Option<usize> {
		None
	}
}

impl AbstractSatSolver for CaDiCaL<'static, 'static> {
	fn clause_count(&self) -> usize {
		self.n_clauses()
	}
	fn label(&self) -> &'static str {
		"CaDiCaL"
	}
}

impl AbstractSatSolver for BatSat {
	fn clause_count(&self) -> usize {
		self.n_clauses()
	}
	fn label(&self) -> &'static str {
		"BatSat"
	}
}

impl AbstractSatSolver for Glucose {
	fn clause_count(&self) -> usize {
		self.n_clauses()
	}
	fn label(&self) -> &'static str {
		"Glucose"
	}
	fn cdcl_constraints(&self) -> Option<usize> {
		Some(self.n_learnts() as usize)
	}
}

#[derive(Debug, Clone)]
pub(super) enum TheoryResult<C: Complex> {
	Maybe(Vec<MaybeData<C>>),
	Conflict(Vec<Clause>),
}

#[derive(Debug)]
pub struct QsatInstance<C: Complex, SatSolver: AbstractSatSolver = Glucose> {
	num_qubits: u16,
	max_depth: u16,
	promise_min_sat: Option<C::Real>,
	constraints: Vec<Constraint<C>>,
	sat: PhantomData<SatSolver>,
}

impl <C: Complex, S: AbstractSatSolver> QsatInstance<C, S> {
	pub fn new(num_qubits: u16, max_depth: u16,) -> Self {
		QsatInstance {
			num_qubits,
			max_depth,
			constraints: Vec::new(),
			promise_min_sat: None,
			sat: PhantomData,
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

	fn learn(&self, assignment: Assignment, num_calls_solver: &mut usize) -> TheoryResult<C> {
		let mut maybes = Vec::new();
		let mut conflicts = Vec::new();
		for c in &self.constraints {
			let result = assignment.find_conflict(c, true, num_calls_solver);
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

pub fn solve<C: Complex, S: AbstractSatSolver>(instance: QsatInstance<C, S>) -> QsatResult<C> {
	solve_diagnostic(instance).result
}

pub fn solve_diagnostic<C: Complex, S: AbstractSatSolver>(
	instance: QsatInstance<C, S>
) -> QsatDiagnosis<C> {
	let mut solver = S::default();
	let mut assignment = instance.assignment();
	if instance.num_vars() > 0 {
		solver.reserve(Var::new(instance.num_vars() - 1))
			.unwrap_or_else(|_| panic!("failed to reserve variables in {}", solver.label()));
	}

	let start = ProcessTime::now();

	let mut num_calls_solver = 0;

	let result = loop {
		let mut temp_assignment = assignment.clone();
		temp_assignment.fill();
		match instance.learn(temp_assignment, &mut num_calls_solver) {
			TheoryResult::Maybe(maybes) => {
				let data = MaybeData::merge(maybes);
				break QsatResult::Maybe(data);
			}
			TheoryResult::Conflict(clauses) => {
				for clause in clauses {
					let sat_clause = clause.to_sat_clause(instance.num_qubits);
					solver.add_clause(sat_clause)
						.unwrap_or_else(|_| panic!("failed to add clause to {}", solver.label()));
				}
			}
		}
		match solver.solve().unwrap_or_else(|_| panic!("{} solve failed", solver.label())) {
			SolverResult::Unsat => {
				break QsatResult::Unsat;
			}
			SolverResult::Sat => {}
			SolverResult::Interrupted => panic!("{} solve was interrupted", solver.label()),
		}
		assignment.update_from_solver(&solver, instance.num_vars());
	};

	let duration = start.elapsed();

	QsatDiagnosis {
		result,
		execution_time: duration,
		num_vars: instance.num_vars() as usize,
		num_calls_solver: num_calls_solver,
		num_clauses: solver.clause_count(),
		num_cdcl_constraints: solver.cdcl_constraints(),
	}
}

#[cfg(test)]
mod test {
	use super::*;
	use rustsat::solvers::{Solve, SolverResult};
	use rustsat::types::{Lit, TernaryVal, Var};
	use rustsat_batsat::BasicSolver as BatSat;
	use rustsat_cadical::CaDiCaL;
	use rustsat_glucose::core::Glucose;

	fn test_cdcl_signs_impl<S: Solve + Default>() {
		let mut solver = S::default();
		let var = Var::new(0);
		let lit = Lit::positive(0);
		let clause: rustsat::types::Clause = vec![lit].into_iter().collect();
		solver.add_clause(clause).expect("failed to add unit clause");
		let result = solver.solve().expect("solve failed");
		assert!(matches!(result, SolverResult::Sat));
		assert_eq!(solver.var_val(var).expect("failed to read value"), TernaryVal::True);
	}

	#[test]
	fn test_cdcl_signs_cadical() {
		test_cdcl_signs_impl::<CaDiCaL<'static, 'static>>();
	}

	#[test]
	fn test_cdcl_signs_batsat() {
		test_cdcl_signs_impl::<BatSat>();
	}

	#[test]
	fn test_cdcl_signs_glucose() {
		test_cdcl_signs_impl::<Glucose>();
	}
}
