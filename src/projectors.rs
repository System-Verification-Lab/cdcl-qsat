use crate::{Complex, Real};

pub fn computational<C: Complex>(element: u16, num_qubits: u16) -> Vec<C> {
	let dim = 1 << num_qubits;
	let mut projector = vec![C::zero(); dim];
	projector[element as usize] = C::raw(C::Real::int(1), C::Real::zero());
	projector
}