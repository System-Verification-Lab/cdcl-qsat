use super::ConvexPolygon;
use crate::{Complex, Real, cross_product};

/// Andrew's convex hull algorithm adapted from
/// https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain

pub(super) fn convex_hull<C: Complex>(points: Vec<C>) -> ConvexPolygon<C> 
{
    let mut points: Vec<C> = points
        .into_iter()
        .flat_map(|point| point.enclosure_points())
        .collect();
    points.sort_by(|a, b| a.re().fine_cmp(&b.re()).then_with(|| a.im().fine_cmp(&b.im())));
    points.dedup_by(|left, right| C::same_point(left, right));
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
            let last = &hull[hull.len() - 1];
            let previous = &hull[hull.len() - 2];
            if C::orientation_sign(previous, last, point) > 0 {
                break;
            }
            hull.pop();
        }
        hull.push(point.clone());
    }
    hull
}

/// Sums and products

pub fn minkowski_sum<C: Complex>(
    mut poly_a: ConvexPolygon<C>,
    mut poly_b: ConvexPolygon<C>,
) -> ConvexPolygon<C> {
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
			let axb_sign = axb.sign();
			if axb_sign > 0 {
				i += 1;
			} else if axb_sign < 0 {
				j += 1;
			} else {
				return ConvexPolygon::new(a_points.iter().flat_map(|a| {
					b_points.iter().map(move |b| a.clone() + b.clone())
				}).collect());
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
    use crate::{cf, numbers::ComplexFloat};

    fn assert_contains_points(
        polygon: &ConvexPolygon<ComplexFloat>,
        expected_points: impl IntoIterator<Item = ComplexFloat>,
    ) {
        for expected in expected_points {
            let translated = ConvexPolygon::new(
                polygon
                    .points()
                    .iter()
                    .map(|point| point.clone() - expected.clone())
                    .collect(),
            );
            assert!(
                translated.contains_zero(),
                "polygon does not enclose expected point {:?}",
                expected.to_float_pair(),
            );
        }
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

    #[test]
    fn convex_hull_keeps_an_exact_turn_hidden_by_interval_rounding() {
        let origin = cf(0.0, 0.0);
        let middle = cf(1.0, 1.0);
        let end = cf(2.0, f64::from_bits(2.0f64.to_bits() + 1));

        // Direct outward-rounded arithmetic cannot certify this very small
        // positive determinant, even though the input coordinates are exact
        // binary64 values and therefore have an exact rational orientation.
        let rounded_turn = cross_product(
            middle.clone() - origin.clone(),
            end.clone() - origin.clone(),
        );
        assert_eq!(rounded_turn.sign(), 0);

        let points = vec![origin, middle, end];
        let hull = convex_hull(points);
        assert_eq!(hull.points().len(), 3);
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
        assert_contains_points(&sum, expected);
    }

    #[test]
    fn minkowski_sum_preserves_zero_when_edge_order_is_uncertain() {
        let left = ConvexPolygon::new(vec![cf(0.0, 0.0), cf(1.0, 1.0)]);
        let right = ConvexPolygon::new(vec![
            cf(-1.0, -1.0),
            cf(0.0, -3.0 * f64::EPSILON),
        ]);

        // (1, 1) + (-1, -1) is exactly zero. The first edges are almost,
        // but not exactly, parallel, so both must not be skipped at once.
        assert!(minkowski_sum(left, right).contains_zero());
    }

}
