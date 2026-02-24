#![allow(unused)]
mod numbers;
mod checker;
mod solver;
mod projectors;
mod utils;
pub use numbers::{Real, Complex, cross_product, f, cf};
pub use solver::{MaybeData};


use numbers::{ComplexFloat};
use projectors::{random_state, random_support};
use solver::{QsatInstance, solve, QsatResult, Glucose};


// Basic usage:
// - Build a QSAT instance with a qubit count and search depth.
// - Add constraints as (support, projector) pairs.
// - Call `solve`` and inspect the result.

type NumberType = ComplexFloat;
type SatSolver = Glucose;

fn main() {
	// Problem parameters.
	let num_qubits = 3;
	let k = 3;
	let num_constraints = 4;
	let search_depth = 8;

	// Create a new QSAT instance with the chosen numeric type.
	let mut instance = QsatInstance::<NumberType, SatSolver>::new(num_qubits, search_depth);

	let seed = 42;
	let max_denominator = 1 << (search_depth);
	let num_amplitudes = 1 << k;

	// Each constraint is passed as a (support, projector) pair.
	// - `support` identifies the qubits involved
	// - `projector` contains the complex conjugate of 2^k amplitudes of the projector state.
	for i in 0..num_constraints {
		let random_support = random_support(num_qubits, k, seed + 7 * i);
		let random_state = random_state::<NumberType>(
			k as u16, num_amplitudes, max_denominator,
			seed + 13 * i, true
		);
		instance.add_constraint(random_support, random_state);
	}

	// Solve and check the outcome.
	let result = solve(instance);
	println!("Result: {:?}", result);
}
