use crate::numbers::{ComplexFloat, Complex, Float};
use crate::projectors::{random_state, random_support};
use crate::solver::{QsatInstance, solve_diagnostic, QsatResult, QsatDiagnosis};

use std::thread;
use std::sync::{Mutex, Arc};

pub fn sample_case<C: Complex>(
    num_qubits: u16, k: u32, num_constraints: usize, promise_gap: C::Real, seed: u64
) -> QsatDiagnosis {

    let mut instance = QsatInstance::<C>::new(num_qubits, promise_gap, k);
    let max_denominator = 8;
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

struct Orchestrator {
    end_qubit: u16,
    num_samples: u32,
    current_qubit: u16,
    current_constraints: usize,
    current_sample: u32
}

impl Orchestrator {
    fn new(start_qubit: u16, end_qubit: u16, num_samples: u32) -> Self {
        Orchestrator {
            end_qubit,
            num_samples,
            current_qubit: start_qubit,
            current_constraints: start_qubit as usize + 1,
            current_sample: 0
        }
    }
    fn next(&mut self) -> Option<(u16, usize, u32)> {
        if self.current_qubit > self.end_qubit {
            return None;
        }
        let result = (self.current_qubit, self.current_constraints, self.current_sample);
        self.current_sample += 1;
        if self.current_sample >= self.num_samples {
            self.current_sample = 0;
            self.current_constraints += 1;
            if self.current_constraints > (self.current_qubit + 1) as usize {
                self.current_qubit += 1;
                self.current_constraints = self.current_qubit as usize + 1;
            }
        }
        Some(result)
    }
}


pub fn benchmark() {
    let k = 3;
    let pgap = 0.125/512.0;

    println!("num_qubits, num_constraints, unprodsat, time, num_calls_solver, num_clauses, num_cdcl_constraints, num_vars, sample_id");

    let orchestrator = Arc::new(Mutex::new(Orchestrator::new(3, 9, 128)));
    let threads = 64;

    let handles: Vec<_> = (0..threads)
        .map(|_| {
            let orchestrator = Arc::clone(&orchestrator);
            thread::spawn(move || {
                loop {
                    let (num_qubits, num_constraints, sample_id) = {
                        let mut orchestrator = orchestrator.lock().unwrap();
                        let Some(next) = orchestrator.next() else { break; };
                        next
                    };
                    benchmark_inner(num_qubits, num_constraints, k, pgap, sample_id)
                }
            })
        })
        .collect();

    for handle in handles {
        handle.join().unwrap();
    }
}

fn benchmark_inner(num_qubits: u16, num_constraints: usize, k: u32, promise_gap: f64, sample_id: u32) {
    let seed = 42
        + (num_qubits as u64) * 32
        + (num_constraints as u64) * 75
        + (sample_id as u64 + 64);

    let diagnosis = sample_case::<ComplexFloat>(
        num_qubits as u16, k, num_constraints, Float::singleton(promise_gap),seed 
    );
    let nma = if let QsatResult::NeedMoreAccuracy = diagnosis.result {
        "!NMA!"
    } else { "" };
    let is_unsat = matches!(diagnosis.result, QsatResult::Unsat);
    let num_vars = diagnosis.num_vars;
    let time = diagnosis.execution_time;
    println!("{}{}, {}, {}, {}, {}, {}, {}, {}, {}",
        nma, num_qubits, num_constraints, is_unsat as u16, time.as_secs_f64(), diagnosis.num_checker_calls, diagnosis.num_clauses, diagnosis.num_cdcl_constraints.unwrap_or(0), num_vars, sample_id);

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn main_example_does_not_report_false_unsat() {
        let diagnosis = sample_case::<ComplexFloat>(
            3,
            3,
            1,
            Float::singleton(0.125 / 64.0),
            231,
        );
        assert!(matches!(diagnosis.result, QsatResult::Sat));
    }
}
