use super::{Assignment, Variable, Clause};
use rustsat::solvers::Solve;
use rustsat::types::{Clause as SatClause, Lit, TernaryVal, Var};

fn to_bool(value: TernaryVal) -> Option<bool> {
    match value {
        TernaryVal::DontCare => None,
        TernaryVal::True => Some(true),
        TernaryVal::False => Some(false),
    }
}

fn to_lit(var: Variable, value: bool, num_qubits: u16) -> Lit {
    if value {
        Lit::positive(var.index(num_qubits))
    } else {
        Lit::negative(var.index(num_qubits))
    }
}

impl Assignment {
    pub(super) fn update_from_solver<S: Solve>(&mut self, solver: &S, num_vars: u32) {
        if num_vars == 0 {
            return;
        }
        let nqubits = self.num_qubits();
        let model = solver.solution(Var::new(num_vars - 1))
            .expect("failed to read solution from SAT solver");
        for var_idx in 0..num_vars as usize {
            let val = model.var_value(Var::new(var_idx as u32));
            let var = Variable::from_index(var_idx, nqubits);
            self.assign(var, to_bool(val))
        }
    }
}

impl Clause {
    pub(super) fn to_sat_clause(self, nqubits: u16) -> SatClause {
        self.0.iter().map(|(var, val)| to_lit(*var, *val, nqubits)).collect()
    }
}
