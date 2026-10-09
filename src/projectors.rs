use crate::{Complex, Real, ComplexFloat, Float, cf};

pub fn computational<C: Complex>(element: u16, num_qubits: u16) -> Vec<C> {
    let dim = 1 << num_qubits;
    let mut projector = vec![C::zero(); dim];
    projector[element as usize] = C::raw(C::Real::int(1), C::Real::zero());
    projector
}

pub fn combine<C: Complex>(a: &Vec<C>, b: &Vec<C>, weight: (u32, u32), phase: (u32, u32)) -> Vec<C> {
    let dim = a.len();
    let weight = (weight.0 as usize, weight.1 as usize);
    assert!(a.len() == b.len());
    let mut combined = Vec::with_capacity(dim);
    let weight_a = C::from_real(C::Real::sin(weight.0, weight.1));
    let weight_b = C::from_real(C::Real::cos(weight.0, weight.1));
    let phase = C::raw(
        C::Real::cos(phase.0 as usize, phase.1 as usize),
        C::Real::sin(phase.0 as usize, phase.1 as usize),
    );
    for i in 0..dim {
        let next = a[i].clone() * weight_a.clone() * phase.clone()
            + b[i].clone() * weight_b.clone();
        combined.push(next);
    }
    combined
}

// Random state generation.
use rand::{SeedableRng, RngExt};
use rand::seq::SliceRandom;
use rand::{rngs::StdRng};


pub fn random_state<C: Complex>(
    num_qubits: u16, num_amplitudes: u32, max_denominator: u32, seed: u64, always_max_den: bool)
 -> Vec<C> {
    assert!(max_denominator > 0, "max_denominator must be positive");
    // Choose `num_amplitudes` random basis states (0..2^numqubits) and store them in `states`.
    let dim = 1 << num_qubits;
    assert!(num_amplitudes <= dim as u32);
    

    let mut rng = StdRng::seed_from_u64(seed);
    let mut indices: Vec<usize> = (0..dim).collect();
    indices.shuffle(&mut rng);
    indices.truncate(num_amplitudes as usize);
    let mut states = Vec::with_capacity(indices.len());
    for &index in indices.iter() {
        states.push(computational::<C>(index as u16, num_qubits));
    }

    // Choose `num_amplitudes - 1` random fractions (num, den) for angles (num < den).
    let mut angles = Vec::with_capacity(num_amplitudes.saturating_sub(1) as usize);
    let mut phase_angles = Vec::with_capacity(num_amplitudes.saturating_sub(1) as usize);
    for _ in 0..num_amplitudes.saturating_sub(1) {
        let den = if always_max_den { max_denominator } else { rng.random_range(1..=max_denominator) };
        let num = rng.random_range(0..den);
        angles.push((num, den));

        let den = if always_max_den { max_denominator } else { rng.random_range(1..=max_denominator) };
        let num = rng.random_range(0..den);
        phase_angles.push((num, den));
    }

    // Combine
    if states.is_empty() {
        return vec![C::zero(); dim];
    }
    let mut state = states[0].clone();
    for (next_state, (angle, phase_angle)) in states.iter().skip(1).zip(angles.iter().zip(phase_angles.iter())) {
        state = combine(&state, next_state, *angle, *phase_angle);
    }

    state
}

#[allow(unused)]
pub fn random_float_state(num_qubits: u16, seed: u64) -> Vec<ComplexFloat> {
    let n = 1 << num_qubits;
    let mut rng = StdRng::seed_from_u64(seed);
    let mut v = Vec::with_capacity(n);

    loop {
        for _ in 0..n {
            v.push(cf(
                rng.random_range(-1.0..1.0),
                rng.random_range(-1.0..1.0),
            ));
        }

        let norm_sq: Float = v.iter()
            .map(|z| z.abs2())
            .reduce(|a, b| a + b)
            .unwrap_or_else(|| Float::zero());

        if norm_sq > Float::zero() {
            let norm = norm_sq.sqrt();

            for z in &mut v {
                *z = z.clone() / ComplexFloat::from_real(norm.clone());
            }

            return v;
        }
    }

}

pub fn random_support(num_qubits: u16, support_size: u32, seed: u64)
 -> Vec<u16> {
    // Return a vector of unique integers between 0 and num_qubits.
    let mut rng = StdRng::seed_from_u64(seed);
    let mut indices: Vec<u16> = (0..num_qubits).collect();
    indices.shuffle(&mut rng);
    indices.truncate(support_size as usize);
    indices
 }
