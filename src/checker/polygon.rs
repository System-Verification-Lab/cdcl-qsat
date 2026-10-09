use super::minkowski::convex_hull;
use crate::{Complex, Real};
use std::cmp::Ordering;

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
            if vec[i].re().fine_cmp(&min_x) == Ordering::Less ||
            (
                vec[i].re().fine_cmp(&min_x) == Ordering::Equal
                 && vec[i].im().fine_cmp(&min_y) == Ordering::Less
            ) {
                min_x = vec[i].re();
                min_y = vec[i].im();
                min_index = i;
            }
        }
        self.0.rotate_left(min_index);
    }

    pub fn normalise(&mut self) {
        let n = self.0.len();

        if n <= 2 {
            self.sort();
            return;
        }

        let reverse = (0..n).find_map(|i| {
            let turn = C::orientation_sign(
                &self.0[i],
                &self.0[(i + 1) % n],
                &self.0[(i + 2) % n],
            );

            match turn {
                -1 => Some(true),  // Clockwise: reverse it.
                1 => Some(false),  // Counter-clockwise: keep it.
                _ => None,         // Collinear: examine another turn.
            }
        });

        if reverse == Some(true) {
            self.0.reverse();
        }

        // Re-establish the canonical starting vertex.
        self.sort();
    }

    pub fn contains_zero(&self) -> bool {
        let vertices = &self.0;
        if vertices.len() == 1 {
            return point_maybe_zero(&vertices[0]);
        } else if vertices.len() == 2 {
            return segment_maybe_contains_zero(&vertices[0], &vertices[1]);
        } else {
            let n = vertices.len();
            let mut sign = 0;
            let mut all_collinear_with_zero = true;
            for i in 0..n {
                let new_sign = C::cross_sign(&vertices[i], &vertices[(i + 1) % n]);
                all_collinear_with_zero &= new_sign == 0;
                if sign == 0 {
                    sign = new_sign;
                } else if new_sign != 0 && new_sign != sign {
                    return false;
                }
            }
            if all_collinear_with_zero {
                return (0..n)
                    .any(|i| segment_maybe_contains_zero(&vertices[i], &vertices[(i + 1) % n]));
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
            } else if !(mod2 < max_mod2) {
                let difference = (max_mod2.clone() - mod2.clone()).abs();
                max_mod2 = (max_mod2 + mod2 + difference) / C::Real::int(2);
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

fn point_maybe_zero<C: Complex>(point: &C) -> bool {
    point.re().sign() == 0 && point.im().sign() == 0
}

fn segment_maybe_contains_zero<C: Complex>(left: &C, right: &C) -> bool {
    point_maybe_zero(left)
        || point_maybe_zero(right)
        || (C::cross_sign(left, right) == 0 && C::dot_sign(left, right) <= 0)
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

    #[test]
    fn max_mod2_encloses_larger_norm_with_overlapping_bounds() {
        let x = 1.0_f64.next_up();
        let y = 1e-8;
        let polygon = ConvexPolygon::raw(vec![cf(1.0, 0.0), cf(x, y)]);
        let (_, upper) = polygon.max_mod2().bounds();
        let exact_x = rug::Rational::from_f64(x).unwrap();
        let exact_y = rug::Rational::from_f64(y).unwrap();
        let exact_norm = exact_x.clone() * exact_x + exact_y.clone() * exact_y;

        assert!(rug::Rational::from_f64(upper).unwrap() >= exact_norm);
    }
}
