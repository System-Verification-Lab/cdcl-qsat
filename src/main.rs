#![allow(unused)]
mod numbers;
mod checker;
mod solver;
mod projectors;
mod utils;
pub use numbers::{Real, Complex, cross_product, f, cf};
pub use solver::{MaybeData};


use numbers::{ComplexFloat, ExactComplex};
use projectors::{random_state, random_support};
use solver::{QsatInstance, solve, QsatResult};

type NumberType = ComplexFloat;

fn main() {
	let num_qubits = 3;
	let k = 3;
	let num_constraints = 4;
	let search_depth = 8;

	let mut instance = QsatInstance::<ComplexFloat>::new(num_qubits, search_depth);

	let seed = 42;
	let max_denominator = 1 << (search_depth);
	let num_amplitudes = 1 << k;

	for i in 0..num_constraints {
		let random_support = random_support(num_qubits, k, seed + 7 * i);
		let random_state = random_state::<NumberType>(
			k as u16, num_amplitudes, max_denominator,
			seed + 13 * i, true
		);
		instance.add_constraint(random_support, random_state);
	}

	let result = solve(instance);
	println!("Result: {:?}", result);
}

