use crate::numbers::{ComplexFloat, ExactComplex, Complex};
use crate::projectors::{random_state, random_support};
use crate::solver::{QsatInstance, solve_diagnostic, QsatResult, QsatDiagnosis};

use std::time::Duration;
use std::io::{stdout, Write};

static NUM_BENCH_SAMPLES: usize = 13;

pub fn sample_case<C: Complex>(
	num_qubits: u16, k: usize, num_constraints: usize, depth: u16, seed: u64
) -> QsatDiagnosis<C> {

	let mut instance = QsatInstance::<C>::new(num_qubits, depth);
	let max_denominator = 1 << (depth);
	let num_amplitudes = 1 << k;
	for i in 0..num_constraints {
		let seed1 = 2026*seed + 100 + i as u64;
		let seed2 = 1998 *seed + 100 + i as u64;
		instance.add_constraint(
			random_support(num_qubits, k as u32 , seed1),
			random_state::<C>(
				k as u16, num_amplitudes, max_denominator, seed2,true 
			),
		);
	}
	solve_diagnostic(instance)
}



pub fn benchmark() {
	let k = 3;
	let depth = 8;

	println!("num_qubits, num_constraints, unprodsat, time, num_calls_solver, num_clauses, num_vars");

	for num_qubits in 6..11 {
		let mut results_tmp = Vec::new();
		for num_constraints in 1..(num_qubits+2) {
			let mut unsat_count = 0;
			for delta in 0..NUM_BENCH_SAMPLES {
				let diagnosis = sample_case::<ComplexFloat>(
					num_qubits as u16, k, num_constraints, depth, delta as u64
				);
				let is_unsat = matches!(diagnosis.result, QsatResult::Unsat);
				let num_vars = diagnosis.num_vars;
				let time = diagnosis.execution_time;
				if is_unsat {
					unsat_count += 1;
				}
				println!("{}, {}, {}, {}, {}, {}, {}",
					num_qubits, num_constraints, is_unsat as u16, time.as_secs_f64(), diagnosis.num_calls_solver, diagnosis.num_clauses, num_vars);
				
			}
			results_tmp.push(unsat_count);
		}
	}
}