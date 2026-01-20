use super::{Assignment, Variable, Clause};
use batsat::{BasicSolver, Lit, SolverInterface, Var, lbool};

fn to_bool(value: lbool) -> Option<bool> {
	if value == lbool::UNDEF {
		return None;
	} else if value == lbool::TRUE {
		return Some(true);
	} else {
		return Some(false);
	}
}

fn to_lit(var: Variable, value: bool, num_qubits: u16) -> Lit {
	Lit::new(Var::unsafe_from_idx(var.index(num_qubits)), value)
}

impl Assignment {
	pub(super) fn update_from_solver(&mut self, solver: &BasicSolver) {
		//self.reset();
		let nqubits = self.num_qubits();
		for (var_idx, val) in solver.get_model().iter().enumerate(){
			let var = Variable::from_index(var_idx, nqubits);
			self.assign(var, to_bool(*val));
		}
	}
}

impl Clause {
	pub(super) fn to_sat_clause(self, nqubits: u16) -> Vec<Lit> {
		self.0.iter().map(|(var, val)| to_lit(*var, *val, nqubits)).collect()
	}
}