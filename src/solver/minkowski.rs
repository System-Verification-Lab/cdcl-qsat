use crate::{Complex, Real, cross_product};
use super::ConvexPolygon;

/// Andrew's convex hull algorithm adapted from
/// https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain

pub fn convex_hull<C: Complex>(mut points: Vec<C>) -> ConvexPolygon<C> 
{
	points.sort_by(|a, b| a.re().cmp(&b.re()).then_with(|| a.im().cmp(&b.im())));
	points.dedup();
	if points.len() < 3 {
		return ConvexPolygon::raw(points);
	}
	
	let mut lower_hull = build_hull(points.iter());
	let mut upper_hull = build_hull(points.iter().rev());

	lower_hull.pop();
	upper_hull.pop();
	lower_hull.append(&mut upper_hull);
	ConvexPolygon::raw(lower_hull)
}

fn build_hull<'a, C: Complex + 'a>(points: impl Iterator<Item = &'a C>) -> Vec<C > {
	let mut hull = Vec::<C>::new();
	for point in points {
		while hull.len() >= 2 {
			let p1 = hull[hull.len() - 1].clone();
			let p2 = hull[hull.len() - 2].clone();
			if cross_product(p1 - p2.clone(), point.clone() - p2) > C::Real::zero() {
				break;
			}
			hull.pop();
		}
		hull.push(point.clone());
	}
	hull
}

#[cfg(test)]
mod test {
	use super::convex_hull;
	use crate::arithmetic::{ComplexFloat, Float};
	use crate::Complex;

	fn f(value: f64) -> Float {
		Float::from(value)
	}

	fn cf(re: f64, im: f64) -> ComplexFloat {
		ComplexFloat::raw(f(re), f(im))
	}

	#[test]
	fn convex_hull_dedups_and_sorts_short_input() {
		let points = vec![cf(1.0, 0.0), cf(0.0, 0.0), cf(1.0, 0.0)];
		let hull = convex_hull(points);
		let expected = vec![cf(0.0, 0.0), cf(1.0, 0.0)];
		assert_eq!(hull.points(), &expected);
	}

	#[test]
	fn convex_hull_ignores_interior_points() {
		let points = vec![
			cf(0.0, 0.0),
			cf(0.0, 1.0),
			cf(0.5, 0.5),
			cf(1.0, 0.0),
			cf(1.0, 1.0),
		];
		let hull = convex_hull(points);
		let expected = vec![cf(0.0, 0.0), cf(1.0, 0.0), cf(1.0, 1.0), cf(0.0, 1.0)];
		assert_eq!(hull.points(), &expected);
	}

	#[test]
	fn convex_hull_collapses_collinear_points() {
		let points = vec![cf(0.0, 0.0), cf(1.0, 0.0), cf(2.0, 0.0)];
		let hull = convex_hull(points);
		let expected = vec![cf(0.0, 0.0), cf(2.0, 0.0)];
		assert_eq!(hull.points(), &expected);
	}
}
