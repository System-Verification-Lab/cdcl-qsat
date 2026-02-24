mod semantics;
mod instance;
mod translation;

use semantics::{Clause, ConstraintLearning, Assignment, Variable, Constraint};
pub use instance::{
	MaybeData, QsatDiagnosis, QsatInstance, QsatResult, AbstractSatSolver, solve, solve_diagnostic,
};

pub use rustsat_batsat::BasicSolver as BatSat;
pub use rustsat_cadical::CaDiCaL as CaDiCaL_core;
pub use rustsat_glucose::core::Glucose as Glucose;

pub type CaDiCal = CaDiCaL_core<'static, 'static>;