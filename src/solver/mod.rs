mod minkowski;
mod annulus;
mod zerotest;

use crate::Complex;

struct ConvexPolygon<C: Complex>(Vec<C>);

impl<C: Complex> ConvexPolygon<C> {
	fn raw(points: Vec<C>) -> Self {
		ConvexPolygon(points)
	}
	fn points(&self) -> &Vec<C> {
		&self.0
	}
}
