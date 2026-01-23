mod numbers;
mod checker;
mod solver;
mod projectors;

pub use numbers::{Real, Complex, cross_product, f, cf};
pub use solver::{MaybeData};

use numbers::{ComplexFloat};
use projectors::{computational, combine, random_state, random_support};
use solver::{QsatInstance, solve, QsatResult};

use std::time::Duration;
use cpu_time::ProcessTime;
use std::io::{stdout, Write};

fn sample_case(num_qubits: u16, k: usize, num_constraints: usize, seed: u64)
-> QsatResult<ComplexFloat> {

    let mut instance = QsatInstance::<ComplexFloat>::new(num_qubits,8);
    let num_amplitudes = 1 << k;
    for i in 0..num_constraints {
        let seed1 = 2026*seed + 100 + i as u64;
        let seed2 = 1998 *seed + 100 + i as u64;
        instance.add_constraint(
            random_support(num_qubits, k as u32 , seed1),
            random_state::<ComplexFloat>(
                k as u16, num_amplitudes, 128, seed2
            ),
        );
    }
    solve(instance)
}

fn bench_sample(num_qubits: u16, k: usize, num_constraints: usize, seed: u64) -> (Duration, bool) {
    let start = ProcessTime::now();
    let result = sample_case(num_qubits, k, num_constraints, seed);
    let duration = start.elapsed();
    let result = matches!(result, QsatResult::Unsat);
    (duration, result)
}

static NUM_BENCH_SAMPLES: usize = 10;

fn run_benchmark() -> (Vec<Duration>, Vec<Vec<u32>>) {
    let mut times = Vec::new();
    let mut results = Vec::new();
    for num_qubits in 3..12 {
        println!(">> Benchmarking num_qubits = {}", num_qubits);
        let mut results_tmp = Vec::new();
        for num_constraints in 1..num_qubits {
            println!("num_constraints = {}", num_constraints);
            let mut unsat_count = 0;
            for delta in 0..NUM_BENCH_SAMPLES {
                let (time, is_unsat) = bench_sample(num_qubits as u16, 3, num_constraints, delta as u64);
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

fn benchmark() {
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

fn main() {

    let num_qubits = 6;
    let num_constraints = 5;
    let k = 3;

    for delta in 0..1 {
        let mut instance = QsatInstance::<ComplexFloat>::new(num_qubits,8);

        for i in 0..num_constraints {
            instance.add_constraint(
                random_support(num_qubits, k , 2026 + delta*100 + i as u64),
                random_state::<ComplexFloat>(
                    k as u16, 8, 64, 1998 + delta*100 + i as u64
                ),
            );
        }

        let result = solve(instance);
        println!("Result: {:?}", result);
    }
}

