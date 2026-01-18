use crate::{Complex, Real, cross_product};
use super::ConvexPolygon;

/// Andrew's convex hull algorithm adapted from
/// https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain

pub(super) fn convex_hull<C: Complex>(mut points: Vec<C>) -> ConvexPolygon<C> 
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

/// Sums and products

pub fn convex_minkowski_product<C: Complex>(poly_a: ConvexPolygon<C>, poly_b: ConvexPolygon<C>)
-> ConvexPolygon<C> {
	let mut points = Vec::new();
	for a in poly_a.points() {
		for b in poly_b.points() {
			points.push(a.clone() * b.clone());
		}
	}
	convex_hull(points)
}

pub fn minkowski_sum<C: Complex>(mut poly_a: ConvexPolygon<C>, mut poly_b: ConvexPolygon<C>)
-> ConvexPolygon<C> {
	poly_a.sort(); poly_a.normalise();
	poly_b.sort(); poly_b.normalise();
	let a_points = poly_a.points();
	let b_points = poly_b.points();
	let n = a_points.len();
	let m = b_points.len();
	let mut i = 0;
	let mut j = 0;
	let mut result_points = Vec::with_capacity(n + m);
	while i < n || j < m {
		result_points.push(a_points[i % n].clone() + b_points[j % m].clone());
		if i < n && j < m {
			let a_edge = a_points[(i + 1) % n].clone() - a_points[i % n].clone();
			let b_edge = b_points[(j + 1) % m].clone() - b_points[j % m].clone();
			let axb = cross_product(a_edge.clone(), b_edge.clone());
			if axb > C::Real::zero() {
				i += 1;
			} else if axb < C::Real::zero() {
				j += 1;
			} else {
				i += 1;
				j += 1;
			}
		} else if i < n {
			i += 1;
		} else {
			j += 1;
		}
	}

	ConvexPolygon::new(result_points)
}

#[cfg(test)]
mod test {
	use super::*;
	use crate::{cf};

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

	#[test]
	fn minkowski_sum_merges_parallel_edges() {
		let square = ConvexPolygon::raw(vec![
			cf(0.0, 0.0),
			cf(1.0, 0.0),
			cf(1.0, 1.0),
			cf(0.0, 1.0),
		]);
		let rectangle = ConvexPolygon::raw(vec![
			cf(0.0, 0.0),
			cf(2.0, 0.0),
			cf(2.0, 1.0),
			cf(0.0, 1.0),
		]);
		let sum = minkowski_sum(square, rectangle);
		let expected = vec![cf(0.0, 0.0), cf(3.0, 0.0), cf(3.0, 2.0), cf(0.0, 2.0)];
		assert_eq!(sum.points(), &expected);
	}

	#[test]
	fn convex_minkowski_product_hulls_pairwise_products() {
		let triangle = ConvexPolygon::raw(vec![
			cf(0.0, 0.0),
			cf(1.0, 0.0),
			cf(0.0, 1.0),
		]);
		let hull = convex_minkowski_product(triangle.clone(), triangle);
		let expected = vec![cf(-1.0, 0.0), cf(1.0, 0.0), cf(0.0, 1.0)];
		assert_eq!(hull.points(), &expected);
	}

}
