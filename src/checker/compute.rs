use crate::numbers::{Complex, Real};
use super::ConvexPolygon;
use super::minkowski::{convex_minkowski_product, minkowski_sum};
use crate::MaybeData;

use std::fmt::Debug;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AngleFraction {
	numerator: usize,
	denominator: usize,
}

impl AngleFraction {
	pub fn new(numerator: usize, denominator: usize) -> Self {
		AngleFraction { numerator, denominator }
	}
	pub fn smaller_one(&self) -> bool {
		self.numerator < self.denominator
	}
	pub fn scale(self, factor: usize) -> Self {
		AngleFraction {
			numerator: self.numerator * factor,
			denominator: self.denominator * factor,
		}
	}
}

pub(super) fn basic_polygons<C: Complex>(theta: AngleFraction, phi: AngleFraction)
-> [ConvexPolygon<C>; 2]
{
	if !theta.smaller_one() || !phi.smaller_one() {
		panic!("Angles must be less than 1");
	}

	// Theta sines/cosines.
	let cos_t0 = C::Real::cos(theta.numerator, theta.denominator);
	let cos_t1 = C::Real::cos(theta.numerator + 1, theta.denominator);
	let sin_t0 = C::Real::sin(theta.numerator, theta.denominator);
	let sin_t1 = C::Real::sin(theta.numerator + 1, theta.denominator);

	// Phi sines/cosines, doubling the resolution to find the mid angle.

	let phi = phi.scale(2);

	let cos_p0 = C::Real::cos(phi.numerator, phi.denominator);
	let cos_p1 = C::Real::cos(phi.numerator + 2, phi.denominator);
	let sin_p0 = C::Real::sin(phi.numerator, phi.denominator);
	let sin_p1 = C::Real::sin(phi.numerator + 2, phi.denominator);
	let cos_pmid = C::Real::cos(phi.numerator + 1, phi.denominator);
	let sin_pmid = C::Real::sin(phi.numerator + 1, phi.denominator);

	let cos_pdelta = C::Real::cos(1, phi.denominator);

	// Build polygons.
	let a = C::from_real(sin_t1.clone()) *C::raw(cos_pmid.clone(), sin_pmid.clone());
	let a = a / C::from_real(cos_pdelta.clone());

	let b = C::from_real(sin_t1.clone()) * C::raw(cos_p0.clone(), sin_p0.clone());
	let c = C::from_real(sin_t0.clone()) * C::raw(cos_p0.clone(), sin_p0.clone());
	let d = C::from_real(sin_t0.clone()) * C::raw(cos_p1.clone(), sin_p1.clone());
	let e = C::from_real(sin_t1.clone()) * C::raw(cos_p1.clone(), sin_p1.clone());

	let interval = ConvexPolygon::<C>::new(vec![C::from_real(cos_t0), C::from_real(cos_t1)]);
	let annulus = ConvexPolygon::<C>::new(vec![a, b, c, d, e]);

	[interval, annulus]

}



#[derive(Debug, Clone)]
pub enum ConstraintResult<C: Complex> {
	Unsat,
	Maybe(MaybeData<C>),
}

impl<C: Complex> ConstraintResult<C> {
	pub fn is_unsat(&self) -> bool {
		matches!(self, ConstraintResult::Unsat)
	}
	pub fn is_maybe(&self) -> bool {
		matches!(self, ConstraintResult::Maybe(_))
	}
}


pub fn poly_check_constraint<C: Complex>(
	projector: &Vec<C>, angle_limits: Vec<(AngleFraction, AngleFraction)>
) -> (ConstraintResult<C>, ConvexPolygon<C>) {

	let n = angle_limits.len();
	let space_dim = projector.len();
	assert_eq!(space_dim, 2usize.pow(n as u32));

	let polygons: Vec<[ConvexPolygon<C>; 2]> = angle_limits.into_iter()
		.map(|(theta, phi)| basic_polygons::<C>(theta, phi))
		.collect();

	let mut summands = Vec::<ConvexPolygon<C>>::with_capacity(space_dim);
	for (i, coeff) in projector.into_iter().enumerate() {
		let poly = polygons
			.iter().enumerate()
			.map(|(pos, polygons)| polygons[(i >> pos) & 1].clone())
			.reduce(|a, b| convex_minkowski_product(a, b))
			.unwrap();
		summands.push(poly.scale(coeff.clone()));
	}
	let final_cover = summands.into_iter().reduce(|a, b| minkowski_sum(a, b)).unwrap();

	if !final_cover.contains_zero() {
		(ConstraintResult::Unsat, final_cover)
	} else {
		let area = final_cover.area();
		let max_mod2 = final_cover.max_mod2();
		(ConstraintResult::Maybe(MaybeData { area, max_mod2,}), final_cover)
	}
}

pub fn check_constraint<C: Complex>(
	projector: &Vec<C>, angle_limits: Vec<(AngleFraction, AngleFraction)>
) -> ConstraintResult<C> {
	let (result, _poly) = poly_check_constraint::<C>(projector, angle_limits);
	result
}


#[cfg(test)]
mod tests {
	use super::*;
	use crate::cf;

	#[test]
	fn test_constraint_checker() {
		let mut projector = vec![cf(0.0,0.0); 8];
		projector[0] = cf(1./2f64.sqrt(), 0.0);
		projector[7] = cf(-1./2f64.sqrt(), 0.0);

		let resolution = 10;
		let sat_angles = vec![
			(AngleFraction::new(0, resolution), AngleFraction::new(0,resolution)),
			(AngleFraction::new(resolution/4, resolution), AngleFraction::new(0,resolution)),
			(AngleFraction::new(0, resolution), AngleFraction::new(0,resolution)),
		];

		let unsat_angles = vec![
			(AngleFraction::new(0, resolution), AngleFraction::new(0,resolution)),
			(AngleFraction::new(0, resolution), AngleFraction::new(0,resolution)),
			(AngleFraction::new(0, resolution), AngleFraction::new(0,resolution)),
		];

		let (sat, _) = poly_check_constraint(&projector, sat_angles);
		let (unsat, _) = poly_check_constraint(&projector, unsat_angles);
		

		assert!(matches!(sat, ConstraintResult::Maybe { .. }));
		assert!(matches!(unsat, ConstraintResult::Unsat));

		/*
		println!("SAT points");
		sat_poly.print_points();
		println!("");

		println!("\nUNSAT points");
		unsat_poly.print_points();
		println!("");

		*/

	}
}
