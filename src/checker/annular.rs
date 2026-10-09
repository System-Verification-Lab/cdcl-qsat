use crate::{Complex, Real};
use super::ConvexPolygon;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AngleInterval {
    start_numerator: usize,
    span_numerator: usize,
    denominator: usize,
}

impl AngleInterval {
    pub fn new(start_numerator: usize, end_numerator: usize, denominator: usize) -> Self {
        assert!(denominator > 0, "angle denominator must be positive");
        assert!(
            start_numerator < denominator && end_numerator <= denominator,
            "angle numerators must be less than denominator"
        );
        let span_numerator = if end_numerator >= start_numerator {
            end_numerator - start_numerator
        } else {
            denominator - start_numerator + end_numerator
        };
        AngleInterval { start_numerator, span_numerator, denominator }
    }

    pub fn from_start_and_span(
        start_numerator: usize,
        span_numerator: usize,
        denominator: usize,
    ) -> Self {
        assert!(denominator > 0, "angle denominator must be positive");
        assert!(start_numerator < denominator, "start_numerator must be less than denominator");
        assert!(span_numerator <= denominator, "span_numerator cannot exceed one full turn");
        AngleInterval { start_numerator, span_numerator, denominator }
    }

    pub fn full_turn(denominator: usize) -> Self {
        assert!(denominator > 0, "angle denominator must be positive");
        AngleInterval {
            start_numerator: 0,
            span_numerator: denominator,
            denominator,
        }
    }

    pub fn start_numerator(&self) -> usize {
        self.start_numerator
    }

    pub fn end_numerator(&self) -> usize {
        (self.start_numerator + self.span_numerator) % self.denominator
    }

    pub fn span_numerator(&self) -> usize {
        self.span_numerator
    }

    pub fn denominator(&self) -> usize {
        self.denominator
    }

    fn rescale(&self, denominator: usize) -> Self {
        assert!(denominator % self.denominator == 0, "new denominator must be a multiple");
        let factor = denominator / self.denominator;
        AngleInterval {
            start_numerator: self.start_numerator * factor,
            span_numerator: self.span_numerator * factor,
            denominator,
        }
    }

    fn add(&self, other: &Self) -> Self {
        let denominator = lcm(self.denominator, other.denominator);
        let left = self.rescale(denominator);
        let right = other.rescale(denominator);
        let start_numerator = (left.start_numerator + right.start_numerator) % denominator;
        let span_numerator = (left.span_numerator + right.span_numerator).min(denominator);
        AngleInterval {
            start_numerator,
            span_numerator,
            denominator,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AnnularSection<C: Complex> {
    inner_radius: C::Real,
    outer_radius: C::Real,
    angles: AngleInterval,
}

impl<C: Complex> AnnularSection<C> {
    pub fn new(inner_radius: C::Real, outer_radius: C::Real, angles: AngleInterval) -> Self {
        assert!(inner_radius.sign() >= 0, "inner_radius must be non-negative");
        assert!(outer_radius.sign() >= 0, "outer_radius must be non-negative");
        assert!(inner_radius <= outer_radius, "inner_radius must be <= outer_radius");
        AnnularSection {
            inner_radius,
            outer_radius,
            angles,
        }
    }

    pub fn segment(r0: C::Real, r1: C::Real) -> Self {
        Self::from_angle_bounds(r0, r1,
            0, 0, 1)
    }

    pub fn from_angle_bounds(
        inner_radius: C::Real,
        outer_radius: C::Real,
        start_numerator: usize,
        end_numerator: usize,
        denominator: usize,
    ) -> Self {
        Self::new(
            inner_radius,
            outer_radius,
            AngleInterval::new(start_numerator, end_numerator, denominator),
        )
    }

    pub fn inner_radius(&self) -> &C::Real {
        &self.inner_radius
    }

    pub fn outer_radius(&self) -> &C::Real {
        &self.outer_radius
    }

    pub fn angles(&self) -> AngleInterval {
        self.angles
    }

    pub fn contains_zero(&self) -> bool {
        self.inner_radius.sign() == 0
    }

    pub fn enclosing_polygon(
        &self,
    ) -> ConvexPolygon<C> {
        let num = self.angles.span_numerator();
        let den = self.angles.denominator();
        let factor = 4;
        let auto_segments = {
            if num > den / 2 { 16 * factor}
            else if num > den / 4 { 8 * factor}
            else if num > den / 8 { 4 * factor}
            else if num > den / 16 { 2 * factor}
            else { 1 * factor }
        };
        let auto_segments = auto_segments.max(1);
        self.enclosing_polygon_with_segments(auto_segments)
    }

    pub fn enclosing_polygon_with_segments(
        &self,
        segments: usize,
    ) -> ConvexPolygon<C> {

        let start = self.angles.start_numerator();
        let end = start + self.angles.span_numerator();
        let denominator = self.angles.denominator();

        let mut points = Vec::with_capacity(segments + 3);

        let start_dir = C::dir(start, denominator);
        let end_dir = C::dir(end, denominator);
        points.push(C::polar(self.inner_radius.clone(), end_dir.clone()));

        let scaled_denominator = denominator * segments * 2;
        let scaled_start = start * segments * 2;
        let half_delta = self.angles.span_numerator();

        if half_delta != 0 {
            points.push(C::polar(self.inner_radius.clone(), start_dir.clone()));

            if 4 * half_delta >= scaled_denominator {
                println!("Half delta: {}, scaled denominator: {}", half_delta, scaled_denominator);
                panic!("Too few segments to approximate annular section");
            }

            let scaled_radius = self.outer_radius.clone() / C::Real::cos(
                half_delta,
                scaled_denominator,
            );

            for i in 0..(segments+1) {
                let direction = C::dir(
                    scaled_start + i * 2 * half_delta, scaled_denominator
                );
                points.push(C::polar(scaled_radius.clone(), direction));
            }
        } else {
            // If the angle span is zero, we just add the outer radius point in that direction.
            points.push(C::polar(self.outer_radius.clone(), start_dir.clone()));
        }
        ConvexPolygon::new(points)
    }
}

pub fn minkowski_product<C: Complex>(
    left: &AnnularSection<C>,
    right: &AnnularSection<C>,
) -> AnnularSection<C> {
    AnnularSection::new(
        left.inner_radius.clone() * right.inner_radius.clone(),
        left.outer_radius.clone() * right.outer_radius.clone(),
        left.angles.add(&right.angles),
    )
}

pub fn minkowski_product_all<C: Complex>(
    sections: impl IntoIterator<Item = AnnularSection<C>>,
) -> Option<AnnularSection<C>> {
    sections
        .into_iter()
        .reduce(|acc, section| minkowski_product(&acc, &section))
}

fn gcd(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    if a == 0 { 1 } else { a }
}

fn lcm(a: usize, b: usize) -> usize {
    if a == b { a } else {
        (a / gcd(a, b)) * b
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::numbers::{ComplexFloat, f};

    #[test]
    fn test_angle_interval_add() {
        let a = AngleInterval::new(1, 3, 8);
        let b = AngleInterval::new(2, 5, 8);
        let c = a.add(&b);
        assert_eq!(c.start_numerator(), (1 + 2) % 8);
        assert_eq!(c.span_numerator(), (2 + 3).min(8));
    }

    #[test]
    fn test_annular_section_enclosing_polygon() {
        let section = AnnularSection::<ComplexFloat>::from_angle_bounds(
            f(1.0), f(2.0),
            0,  0, 8
        );
        let polygon = section.enclosing_polygon_with_segments(1);
        polygon.print_points();
    }
}
