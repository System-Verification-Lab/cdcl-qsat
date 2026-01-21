mod numbers;
mod checker;
mod solver;
mod projectors;
//mod cdcl;

pub use numbers::{Real, Complex, cross_product, f, cf};
pub use solver::{MaybeData};

use numbers::{ComplexFloat};
use projectors::{computational, combine, random_state, random_support};

use solver::{QsatInstance, solve};


fn main() {

    let num_qubits = 6;
    let num_constraints = 5;
    let k = 3;

    for delta in 0..100 {
        let mut instance = QsatInstance::<ComplexFloat>::new(num_qubits,8);

        for i in 0..num_constraints {
            instance.add_constraint(
                random_support(num_qubits, k , 2026 + delta*100 + i as u64),
                random_state::<ComplexFloat>(
                    k as u16, 8, 32, 1998 + delta*100 + i as u64
                ),
            );
        }

        let result = solve(instance);
        println!("Result: {:?}", result);
    }


}

