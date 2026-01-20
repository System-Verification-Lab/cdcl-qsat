mod numbers;
mod checker;
mod solver;
mod projectors;
//mod cdcl;

pub use numbers::{Real, Complex, cross_product, f, cf};
pub use solver::{MaybeData};

use numbers::{ComplexFloat};
use projectors::computational;

use solver::{QsatInstance, solve};

fn main() {

    
    let nqubits: u16 = 3;
    
    let mut instance = QsatInstance::<ComplexFloat>::new(nqubits,3);
    for qubit in 0..(2u16.pow(nqubits as u32)) {
        instance.add_constraint(
            ((0..nqubits)).collect(), computational(qubit, nqubits)
        );
    }
    println!("Instance: {:#?}", instance);
    let result = solve(instance);
    println!("Result: {:?}", result);


}

