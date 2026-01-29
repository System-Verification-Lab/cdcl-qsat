use crate::numbers::{ComplexFloat, ExactComplex, Complex};
use crate::projectors::{random_state, random_support};
use crate::solver::{QsatInstance, solve, QsatResult};

use std::time::Duration;
use cpu_time::ProcessTime;
use std::io::{stdout, Write};

static NUM_BENCH_SAMPLES: usize = 7;

pub fn sample_case<C: Complex>(
	num_qubits: u16, k: usize, num_constraints: usize, depth: u16, seed: u64
) -> QsatResult<C> {

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
	solve(instance)
}

fn bench_sample<C: Complex>(
	num_qubits: u16, k: usize, num_constraints: usize, depth: u16, seed: u64
) -> (Duration, bool) {
	let start = ProcessTime::now();
	let result = sample_case::<C>(num_qubits, k, num_constraints, depth, seed);
	let duration = start.elapsed();
	let result = matches!(result, QsatResult::Unsat);
	(duration, result)
}


fn run_benchmark() -> (Vec<Duration>, Vec<Vec<u32>>) {
	let k = 3;
	let depth = 8;

	let mut times = Vec::new();
	let mut results = Vec::new();
	for num_qubits in 3..5 {
		println!(">> Benchmarking num_qubits = {}", num_qubits);
		let mut results_tmp = Vec::new();
		for num_constraints in 1..(num_qubits+2) {
			println!("num_constraints = {}", num_constraints);
			let mut unsat_count = 0;
			for delta in 0..NUM_BENCH_SAMPLES {
				let (time, is_unsat) = bench_sample::<ComplexFloat>(
					num_qubits as u16, k, num_constraints, depth, delta as u64
				);
				times.push(time);
				if is_unsat {
					unsat_count += 1;
				}
				print!(".");
				stdout().flush().unwrap();
			}
			print!("\n");
			results_tmp.push(unsat_count);
		}
		results.push(results_tmp);
	}
	(times, results)
}

pub fn benchmark() {
	let (times, results) = run_benchmark();
	let time_values: Vec<f64> = times
		.iter()
		.map(|d| d.as_secs_f64())
		.collect(); 
	let avg_time = time_values.iter().sum::<f64>() / time_values.len() as f64;
	let variance = time_values
		.iter()
		.map(|v| {
			let diff = v - avg_time;
			diff * diff
		})
		.sum::<f64>()
		/ time_values.len() as f64;

	let std_time = variance.sqrt();
	let max_time = times.iter().max().cloned().unwrap_or_default();
	println!("Average time: {:?}", avg_time);
	println!("Standard deviation: {:?}", std_time);
	println!("Max time: {:?}", max_time);

	let max_cols = results.iter().map(|r| r.len()).max().unwrap_or(0);
	for row in results.iter() {
		for col in 0..max_cols {
			if col < row.len() {
				print!("{}/{}", row[col], NUM_BENCH_SAMPLES);
			} else {
				print!("");
			}
			print!("&");
		}
		println!("\\\\");
	}
}