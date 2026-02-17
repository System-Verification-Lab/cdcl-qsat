use crate::numbers::{Complex, Real};
use crate::numbers::TrigonometricHelper;
use super::ConvexPolygon;
use super::minkowski::{minkowski_sum};
use super::annular::{AnnularSection, minkowski_product};
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
	pub fn half(self) -> Self {
		AngleFraction {
			numerator: self.numerator,
			denominator: self.denominator * 2
		}
	}
}

pub(super) fn annular_sections<C: Complex>(
	theta: AngleFraction,
	phi: AngleFraction,
)
-> [AnnularSection<C>; 2]
{
	if !theta.smaller_one() || !phi.smaller_one() {
		panic!("Angles must be less than 1");
	}

	let theta = theta.half();

	// Theta sines/cosines.
	let cos_t0 = C::Real::cos(theta.numerator + 1, theta.denominator);
	let cos_t1 = C::Real::cos(theta.numerator, theta.denominator);

	
	let r1 = C::Real::sin(theta.numerator + 1, theta.denominator);
	let r0 = C::Real::sin(theta.numerator, theta.denominator);

	let interval = AnnularSection::<C>::segment(cos_t0, cos_t1);
	let annulus = AnnularSection::<C>::from_angle_bounds(
		r0, r1,
		phi.numerator, phi.numerator + 1, phi.denominator
	);

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
	projector: &Vec<C>,
	angle_limits: Vec<(AngleFraction, AngleFraction)>,
) -> (ConstraintResult<C>, ConvexPolygon<C>) {

	let n = angle_limits.len();
	let space_dim = projector.len();
	assert_eq!(space_dim, 2usize.pow(n as u32));

	let sections: Vec<[AnnularSection<C>; 2]> = angle_limits.into_iter()
		.map(|(theta, phi)| annular_sections::<C>(theta, phi))
		.collect();

	let mut summands = Vec::<ConvexPolygon<C>>::with_capacity(space_dim);
	for (i, coeff) in projector.into_iter().enumerate() {
		let section = sections 
			.iter().enumerate()
			// We use `i >> pos`` because we can equivalently reverse all projectors.
			.map(|(pos, sections)| sections[(i >> pos) & 1].clone())
			.reduce(|a, b| minkowski_product(&a, &b))
			.unwrap();
		summands.push(
			section
				.enclosing_polygon()
				.scale(coeff.clone())
		);
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
	projector: &Vec<C>,
	angle_limits: Vec<(AngleFraction, AngleFraction)>,
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

		let resolution = 16;
		let sat_angles = vec![
			(AngleFraction::new(0, resolution), AngleFraction::new(0,resolution)),
			(AngleFraction::new(resolution/2 - 1, resolution), AngleFraction::new(0,resolution)),
			(AngleFraction::new(0, resolution), AngleFraction::new(0,resolution)),
		];

		let unsat_angles = vec![
			(AngleFraction::new(0, resolution), AngleFraction::new(0,resolution)),
			(AngleFraction::new(0, resolution), AngleFraction::new(0,resolution)),
			(AngleFraction::new(0, resolution), AngleFraction::new(0,resolution)),
		];

		let (sat, _) = poly_check_constraint(&projector, sat_angles,);
		let (unsat, _) = poly_check_constraint(&projector, unsat_angles,);
		

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
