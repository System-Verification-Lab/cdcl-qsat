use core::panic;

use crate::{Complex, Real, cross_product};
use super::minkowski::convex_hull;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConvexPolygon<C: Complex>(Vec<C>);

impl<C: Complex> ConvexPolygon<C> {
	pub fn raw(points: Vec<C>) -> Self {
		if points.is_empty() {
			panic!("ConvexPolygon cannot be empty");
		}
		ConvexPolygon(points)
	}

	pub fn new(points: Vec<C>) -> Self {
		convex_hull(points)
	}

	pub fn points(&self) -> &Vec<C> {
		&self.0
	}

	pub fn sort(&mut self) {
		let vec = &self.0;
		let mut min_index = 0;
		let mut min_x = vec[0].re();
		let mut min_y = vec[0].im();
		for i in 1..vec.len() {
			if vec[i].re() < min_x || (vec[i].re() == min_x && vec[i].im() < min_y) {
				min_x = vec[i].re();
				min_y = vec[i].im();
				min_index = i;
			}
		}
		self.0.rotate_left(min_index);
	}

	pub fn normalise(&mut self) {
		let vertices = &self.0;
		if vertices.len() < 2 {
			return;
		}
		let cross = cross_product(
			vertices[0].clone(),
			vertices[1].clone() - vertices[0].clone());
		if cross < C::Real::zero() {
			self.0.reverse();
		}
	}

	pub fn contains_zero(&self) -> bool {
		let vertices = &self.0;
		if vertices.len() == 1 {
			return vertices[0] == C::zero();
		} else if vertices.len() == 2 {
			let a = vertices[0].clone();
			let b = vertices[1].clone();
			if a == C::zero() || b == C::zero() {
				return true;
			}
			return cross_product(a.clone(), b.clone()) == C::Real::zero()
				&& (a.re() * b.re() < C::Real::zero() || a.im() * b.im() < C::Real::zero())
		} else {
			let n = vertices.len();
			let mut sign = 0;
			for i in 0..n {
				let point = vertices[i].clone();
				let edge = vertices[(i + 1) % n].clone() - point.clone();
				let new_sign = cross_product(point.clone(), edge.clone()).sign();
				if sign == 0 {
					sign = new_sign;
				} else if new_sign != 0 && new_sign != sign {
					return false;
				}
			}
			return true;
		}
	}

	pub fn area(&self) -> C::Real {
		let vertices = &self.0;
		let n = vertices.len();
		let mut area = C::Real::zero();
		for i in 0..n {
			let p1 = vertices[i].clone();
			let p2 = vertices[(i + 1) % n].clone();
			area = area + (p1.re() * p2.im() - p2.re() * p1.im());
		}
		area.abs() / C::Real::int(2)
	}

	pub fn max_mod2(&self) -> C::Real {
		let mut max_mod2 = C::Real::zero();
		for point in &self.0 {
			let mod2 = point.re() * point.re() + point.im() * point.im();
			if mod2 > max_mod2 {
				max_mod2 = mod2;
			}
		}
		max_mod2
	}

	pub fn scale(self, factor: C) -> Self {
		let scaled_points: Vec<C> = self.0.into_iter()
			.map(|point| point * factor.clone())
			.collect();
		ConvexPolygon::raw(scaled_points)
	}
	
	pub fn print_points(&self) {
		print!("[");
		for point in &self.0 {
			let (re, im) = point.to_float_pair();
			print!("[{}, {}],", re, im);
		}
		print!("]");
	}
}

#[cfg(test)]
mod tests {
	use super::ConvexPolygon;
	use crate::numbers::{ComplexFloat, Float};
	use crate::Complex;

	fn f(value: f64) -> Float {
		Float::from(value)
	}

	fn cf(re: f64, im: f64) -> ComplexFloat {
		ComplexFloat::raw(f(re), f(im))
	}

	#[test]
	fn contains_zero_handles_singleton_and_segment() {
		let zero = ConvexPolygon::raw(vec![cf(0.0, 0.0)]);
		let nonzero = ConvexPolygon::raw(vec![cf(1.0, 0.0)]);
		assert!(zero.contains_zero());
		assert!(!nonzero.contains_zero());

		let segment = ConvexPolygon::raw(vec![cf(1.0, 1.0), cf(-2.0, -2.0)]);
		assert!(segment.contains_zero());
	}

	#[test]
	fn contains_zero_checks_polygon_orientation_around_origin() {
		let square = ConvexPolygon::raw(vec![
			cf(1.0, 0.0),
			cf(0.0, 1.0),
			cf(-1.0, 0.0),
			cf(0.0, -1.0),
		]);
		assert!(square.contains_zero());

		let triangle = ConvexPolygon::raw(vec![cf(1.0, 0.0), cf(2.0, 0.0), cf(1.0, 1.0)]);
		assert!(!triangle.contains_zero());
	}
}
