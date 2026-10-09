mod semantics;
mod instance;
mod translation;

use semantics::{Clause, Assignment, Variable, Constraint};
pub use instance::{
    QsatDiagnosis, QsatInstance, QsatResult, AbstractSatSolver, solve, solve_diagnostic,
};

#[allow(unused)]
pub use rustsat_batsat::BasicSolver as BatSat;
pub use rustsat_cadical::CaDiCaL as CaDiCaL_core;

#[allow(unused)]
pub use rustsat_glucose::core::Glucose as Glucose;

#[allow(unused)]
pub type CaDiCal = CaDiCaL_core<'static, 'static>;