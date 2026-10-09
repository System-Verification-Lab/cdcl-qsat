use crate::{Complex, Real};
use super::{Constraint, Assignment, Clause};
use cpu_time::{ThreadTime};
use rustsat::solvers::{Solve, SolveStats, SolverResult};
use rustsat::types::{Var};
use rustsat_batsat::BasicSolver as BatSat;
use rustsat_cadical::CaDiCaL;
use rustsat_glucose::core::Glucose;
use std::time::Duration;
use std::marker::PhantomData;
use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub enum QsatResult {
    Sat,
    Unsat,
    NeedMoreAccuracy,
}

#[derive(Debug, Clone)]
pub struct QsatDiagnosis {
    pub result: QsatResult,
    pub execution_time: Duration,
    pub num_vars: usize,
    pub num_checker_calls: usize,
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
pub(super) enum TheoryResult {
    Sat,
    Conflict(Vec<Clause>),
    NeedMoreAccuracy,
}

#[derive(Debug)]
pub struct QsatInstance<C: Complex, SatSolver: AbstractSatSolver = Glucose> {
    num_qubits: u16,
    max_depth: u16,
    constraints: Vec<Constraint<C>>,
    promise_gap: C::Real,
    support_size: u32,
    sat: PhantomData<SatSolver>,
}

impl <C: Complex, S: AbstractSatSolver> QsatInstance<C, S> {
    pub fn new(num_qubits: u16, promise_gap: C::Real, support_size: u32) -> Self {
        let max_depth = 2.0 * 2f64.sqrt() * PI * support_size as f64 / promise_gap.to_float();
        let max_depth = max_depth.log2().ceil() as u16;
        QsatInstance {
            num_qubits,
            max_depth,
            support_size,
            constraints: Vec::new(),
            sat: PhantomData,
            promise_gap,
        }
    }

    pub fn add_constraint(&mut self, support: Vec<u16>, projector: Vec<C>) {
        let max_qubit = *support.iter().max().unwrap();
        assert!(max_qubit < self.num_qubits, "qubit index out of range");
        assert!(projector.len() == 2usize.pow(support.len() as u32),
            "projector size does not match the number of qubits");
        assert!(support.len() == self.support_size as usize, "support size does not match the expected size");

        self.constraints.push(Constraint::new(support, projector));
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

    fn learn(&self, assignment: Assignment, num_calls_solver: &mut usize) -> TheoryResult {
        let mut found_unsat = vec![];
        for c in 0..self.constraints.len() {
            if self.constraints[c].geometric_unsat(&assignment, num_calls_solver) {
                found_unsat.push(c);
            }
        }
        if found_unsat.is_empty() {
            let mut energy = C::Real::zero();
            for c in &self.constraints {
                energy = energy + c.evaluate(&assignment);
            }
            let energy = energy / C::Real::int(self.constraints.len());
            if energy < self.promise_gap {
                return TheoryResult::Sat;
            } else if !(energy > self.promise_gap.clone() / C::Real::int(2)) {
                return TheoryResult::NeedMoreAccuracy;
            }
        }

        let mut conflicts = Vec::new();

        for c in found_unsat {
            let constraint = &self.constraints[c];
            let clause = assignment.generalise_conflict(constraint, num_calls_solver);
            conflicts.push(clause);
        }
        if conflicts.is_empty() {
            conflicts.push(assignment.full_cell_clause());
        }
        TheoryResult::Conflict(conflicts)
    }

}

pub fn solve<C: Complex, S: AbstractSatSolver>(instance: QsatInstance<C, S>) -> QsatResult {
    solve_diagnostic(instance).result
}

pub fn solve_diagnostic<C: Complex, S: AbstractSatSolver>(
    instance: QsatInstance<C, S>
) -> QsatDiagnosis {
    let mut solver = S::default();
    let mut assignment = instance.assignment();
    if instance.num_vars() > 0 {
        solver.reserve(Var::new(instance.num_vars() - 1))
            .unwrap_or_else(|e| panic!("failed to reserve variables in {}. {e:#?}", solver.label()));
    }

    let start = ThreadTime::now();

    let mut num_calls_solver = 0;

    let result = loop {
        let mut temp_assignment = assignment.clone();
        temp_assignment.fill();
        match instance.learn(temp_assignment, &mut num_calls_solver) {
            TheoryResult::Sat => {
                break QsatResult::Sat;
            }
            TheoryResult::NeedMoreAccuracy => {
                break QsatResult::NeedMoreAccuracy;
            }
            TheoryResult::Conflict(clauses) => {
                for clause in clauses {
                    let sat_clause = clause.to_sat_clause(instance.num_qubits);
                    solver.add_clause(sat_clause)
                        .unwrap_or_else(|e| panic!("failed to add clause to {}. {e:#?}", solver.label()));
                }
            }
        }
        match solver.solve().unwrap_or_else(|e| panic!("{} solve failed. {e:#?}", solver.label())) {
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
        num_checker_calls: num_calls_solver,
        num_clauses: solver.clause_count(),
        num_cdcl_constraints: solver.cdcl_constraints(),
    }
}

#[cfg(test)]
mod test {
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
