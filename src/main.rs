mod numbers;
mod checker;
mod solver;
mod projectors;

pub use numbers::{Real, Complex, cross_product, f, cf};

use solver::{QsatSolver};
use numbers::{ComplexFloat};
use projectors::computational;

use checker::{poly_check_constraint, check_constraint, Fraction};

fn main() {

    
    let nqubits: u16 = 1;
    
    let mut solver = QsatSolver::<ComplexFloat>::new(nqubits,3);
    /*for qubit in 0..(2u16.pow(nqubits as u32)-1) {
        solver = solver.add_constraint(
            ((0..nqubits)).collect(), computational(qubit, nqubits)
        );
    }*/
    //println!("Solver: {:#?}", solver);
    solver = solver.add_constraint(
        vec![0],
        computational::<ComplexFloat>(0, nqubits)
    );
    let result = solver.solve();

    println!("QSAT result: {:?}", result);

    /*
    for th in 0..8 {
        for ph in 0..8 {
            let angle_limits = vec![
                (Fraction::new(th,8), Fraction::new(ph,8)),
            ];

            for state in 0..2 {
                let projector = computational::<ComplexFloat>(state, 1);
                let result = check_constraint(&projector, angle_limits.clone());
                println!("{:?}", result);
            }
            println!("")


        }
    }*/

}

