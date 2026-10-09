use std::cmp::Ordering;
use std::ops::{Add, Div, Mul, Neg, Sub};
use rug::Integer;

use super::{Real, Complex, cross_product};

// Real float interval with outward-rounded endpoints.
#[derive(Clone, Debug)]
pub struct Float {
    lower: f64,
    upper: f64,
}

impl Float {
    pub fn new(lower: impl Into<f64>, upper: impl Into<f64>) -> Self {
        let lower = lower.into();
        let upper = upper.into();
        assert!(
            lower.is_finite() && upper.is_finite(),
            "bounds must be finite"
        );
        assert!(!lower.is_nan() && !upper.is_nan(), "bounds must not be NaN");
        assert!(lower <= upper, "lower bound must be <= upper bound");
        Self::from_bounds(lower, upper)
    }

    pub fn from_bounds(lower: f64, upper: f64) -> Self {
        Self { lower, upper }
    }

    pub fn singleton(value: impl Into<f64>) -> Self {
        let value = value.into();
        Self::new(value, value)
    }

    pub fn bounds(&self) -> (f64, f64) {
        (self.lower, self.upper)
    }

    fn contains(&self, value: f64) -> bool {
        self.lower <= value && value <= self.upper
    }

    fn contains_zero(&self) -> bool {
        self.contains(0.0)
    }

    fn midpoint(&self) -> f64 {
        0.5 * self.lower + 0.5 * self.upper
    }

    pub fn sqrt(&self) -> Self {
        assert!(self.lower >= 0.0, "cannot take square root of negative interval");
        Self::from_bounds(round_down(self.lower.sqrt()), round_up(self.upper.sqrt()))
    }
}

impl From<f64> for Float {
    fn from(value: f64) -> Self {
        Self::singleton(value)
    }
}

impl PartialEq for Float {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Float {}
impl PartialOrd for Float {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Float {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.upper < other.lower {
            Ordering::Less
        }
        else if self.lower > other.upper {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}

impl Add for Float {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self::from_bounds(
            round_down(self.lower + other.lower),
            round_up(self.upper + other.upper),
        )
    }
}

impl Sub for Float {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self::from_bounds(
            round_down(self.lower - other.upper),
            round_up(self.upper - other.lower),
        )
    }
}

impl Mul for Float {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        interval_mul(&self.lower, &self.upper, &other.lower, &other.upper)
    }
}

impl Div for Float {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        assert!(
            !other.contains_zero(),
            "division by interval containing zero"
        );
        interval_div(&self.lower, &self.upper, &other.lower, &other.upper)
    }
}

impl Neg for Float {
    type Output = Self;
    fn neg(self) -> Self {
        Self::from_bounds(-self.upper, -self.lower)
    }
}

impl Real for Float {
    fn zero() -> Self {
        Self::singleton(0.0)
    }

    fn int(from: usize) -> Self {
        let value = from as f64;
        let lower = if (value as u128) > from as u128 {
            value.next_down()
        } else {
            value
        };
        let upper = if (value as u128) < from as u128 {
            value.next_up()
        } else {
            value
        };
        Self::from_bounds(lower, upper)
    }

    fn sign(&self) -> i8 {
        if self.lower > 0.0 {
            1
        } else if self.upper < 0.0 {
            -1
        } else {
            0
        }
    }

    fn abs(self) -> Self {
        if self.lower >= 0.0 {
            self
        } else if self.upper <= 0.0 {
            -self
        } else {
            Self::from_bounds(0.0, (-self.lower).max(self.upper))
        }
    }

    fn cos(numerator: usize, denominator: usize) -> Self {
        trig_bounds(numerator, denominator, TrigFn::Cos)
    }

    fn sin(numerator: usize, denominator: usize) -> Self {
        trig_bounds(numerator, denominator, TrigFn::Sin)
    }

    fn to_float(&self) -> f64 {
        self.midpoint()
    }

    fn fine_cmp(&self, other: &Self) -> Ordering {
        let result = self.cmp(other);
        if result == Ordering::Equal {
            self.midpoint().partial_cmp(&other.midpoint()).unwrap_or(Ordering::Equal)
        } else {
            result
        }
    }
}

// Complex float interval
#[derive(Clone, PartialEq, Debug)]
pub struct ComplexFloat(Float, Float);

impl Add for ComplexFloat {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        ComplexFloat(self.0 + other.0, self.1 + other.1)
    }
}

impl Sub for ComplexFloat {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        ComplexFloat(self.0 - other.0, self.1 - other.1)
    }
}

impl Mul for ComplexFloat {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        ComplexFloat(
            self.0.clone() * other.0.clone() - self.1.clone() * other.1.clone(),
            self.0.clone() * other.1 + self.1.clone() * other.0,
        )
    }
}

impl Div for ComplexFloat {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        let denom = other.0.clone() * other.0.clone() + other.1.clone() * other.1.clone();
        ComplexFloat(
            (self.0.clone() * other.0.clone() + self.1.clone() * other.1.clone()) / denom.clone(),
            (self.1.clone() * other.0 - self.0.clone() * other.1) / denom,
        )
    }
}

impl Neg for ComplexFloat {
    type Output = Self;
    fn neg(self) -> Self {
        ComplexFloat(-self.0, -self.1)
    }
}

impl Complex for ComplexFloat {
    type Real = Float;
    fn re(&self) -> Self::Real {
        self.0.clone()
    }
    fn im(&self) -> Self::Real {
        self.1.clone()
    }
    fn raw(re: Self::Real, im: Self::Real) -> Self {
        ComplexFloat(re, im)
    }
    fn enclosure_points(&self) -> Vec<Self> {
        if self.0.lower == self.0.upper && self.1.lower == self.1.upper {
            return vec![self.clone()];
        }
        let mut points = Vec::with_capacity(4);
        for re in [self.0.lower, self.0.upper] {
            for im in [self.1.lower, self.1.upper] {
                points.push(ComplexFloat(Float::singleton(re), Float::singleton(im)));
            }
        }
        points
    }

    fn same_point(left: &Self, right: &Self) -> bool {
        left.0.lower == right.0.lower
            && left.0.upper == right.0.upper
            && left.1.lower == right.1.lower
            && left.1.upper == right.1.upper
    }

    fn cross_sign(left: &Self, right: &Self) -> i8 {
        match (singleton_coords(left), singleton_coords(right)) {
            (Some((ax, ay)), Some((bx, by))) => sum_of_two_products_sign(
                (ax, by, false),
                (ay, bx, true),
            ),
            _ => cross_product(left.clone(), right.clone()).sign(),
        }
    }

    fn orientation_sign(origin: &Self, left: &Self, right: &Self) -> i8 {
        match (
            singleton_coords(origin),
            singleton_coords(left),
            singleton_coords(right),
        ) {
            (Some((ox, oy)), Some((ax, ay)), Some((bx, by))) => {
                filtered_orientation_sign(ox, oy, ax, ay, bx, by).unwrap_or_else(|| {
                    // Expand the determinant into products of input coordinates,
                    // which are exact dyadic values even when their f64
                    // differences would round or overflow.
                    exact_sum_of_products_sign([
                        (ax, by, false),
                        (ax, oy, true),
                        (ox, by, true),
                        (ay, bx, true),
                        (ay, ox, false),
                        (oy, bx, false),
                    ])
                })
            }
            _ => cross_product(
                left.clone() - origin.clone(),
                right.clone() - origin.clone(),
            )
            .sign(),
        }
    }

    fn dot_sign(left: &Self, right: &Self) -> i8 {
        match (singleton_coords(left), singleton_coords(right)) {
            (Some((ax, ay)), Some((bx, by))) => sum_of_two_products_sign(
                (ax, bx, false),
                (ay, by, false),
            ),
            _ => (left.re() * right.re() + left.im() * right.im()).sign(),
        }
    }
}

fn singleton_coords(value: &ComplexFloat) -> Option<(f64, f64)> {
    if value.0.lower == value.0.upper && value.1.lower == value.1.upper {
        Some((value.0.lower, value.1.lower))
    } else {
        None
    }
}

// A finite binary64 product is an integer of at most 106 bits times a power
// of two. Keeping that representation avoids Rational allocation and GCDs.
#[derive(Clone, Copy)]
struct DyadicProduct {
    negative: bool,
    significand: u128,
    exponent: i32,
}

fn dyadic_product(left: f64, right: f64, negate: bool) -> Option<DyadicProduct> {
    let (left_negative, left_significand, left_exponent) = decompose_finite(left);
    let (right_negative, right_significand, right_exponent) = decompose_finite(right);
    if left_significand == 0 || right_significand == 0 {
        return None;
    }

    Some(DyadicProduct {
        negative: left_negative ^ right_negative ^ negate,
        significand: u128::from(left_significand) * u128::from(right_significand),
        exponent: left_exponent + right_exponent,
    })
}

fn decompose_finite(value: f64) -> (bool, u64, i32) {
    const FRACTION_MASK: u64 = (1_u64 << 52) - 1;

    assert!(value.is_finite(), "geometry coordinates must be finite");
    let bits = value.to_bits();
    let negative = bits >> 63 != 0;
    let raw_exponent = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & FRACTION_MASK;
    if raw_exponent == 0 {
        (negative, fraction, -1074)
    } else {
        (
            negative,
            (1_u64 << 52) | fraction,
            raw_exponent - 1075,
        )
    }
}

fn sum_of_two_products_sign(
    first: (f64, f64, bool),
    second: (f64, f64, bool),
) -> i8 {
    let first = dyadic_product(first.0, first.1, first.2);
    let second = dyadic_product(second.0, second.1, second.2);
    match (first, second) {
        (None, None) => 0,
        (Some(product), None) | (None, Some(product)) => product_sign(product),
        (Some(first), Some(second)) if first.negative == second.negative => {
            product_sign(first)
        }
        (Some(first), Some(second)) => match product_magnitude_cmp(first, second) {
            Ordering::Less => product_sign(second),
            Ordering::Equal => 0,
            Ordering::Greater => product_sign(first),
        },
    }
}

fn product_sign(product: DyadicProduct) -> i8 {
    if product.negative { -1 } else { 1 }
}

fn product_magnitude_cmp(left: DyadicProduct, right: DyadicProduct) -> Ordering {
    let left_bits = u128::BITS - left.significand.leading_zeros();
    let right_bits = u128::BITS - right.significand.leading_zeros();
    let left_top_exponent = left.exponent + left_bits as i32;
    let right_top_exponent = right.exponent + right_bits as i32;

    match left_top_exponent.cmp(&right_top_exponent) {
        // Equal top exponents imply the alignment shift is at most 105 bits
        // and the shifted value still fits in 106 bits.
        Ordering::Equal if left_bits < right_bits => {
            (left.significand << (right_bits - left_bits)).cmp(&right.significand)
        }
        Ordering::Equal if right_bits < left_bits => {
            left.significand.cmp(&(right.significand << (left_bits - right_bits)))
        }
        Ordering::Equal => left.significand.cmp(&right.significand),
        ordering => ordering,
    }
}

fn filtered_orientation_sign(
    ox: f64,
    oy: f64,
    ax: f64,
    ay: f64,
    bx: f64,
    by: f64,
) -> Option<i8> {
    let ax = ax - ox;
    let ay = ay - oy;
    let bx = bx - ox;
    let by = by - oy;
    let left_sign = scalar_sign(ax) * scalar_sign(by);
    let right_sign = scalar_sign(ay) * scalar_sign(bx);

    if left_sign == 0 {
        return Some(-right_sign);
    }
    if right_sign == 0 || left_sign != right_sign {
        return Some(left_sign);
    }

    let left = ax * by;
    let right = ay * bx;
    let determinant = left - right;
    if ![ax, ay, bx, by, left, right, determinant]
        .into_iter()
        .all(f64::is_normal)
    {
        return None;
    }

    // Each product includes two rounded subtractions and one rounded
    // multiplication. Four epsilons times their magnitude sum is a
    // conservative bound on those errors plus the final subtraction.
    let magnitude_sum = (left.abs() + right.abs()).next_up();
    if !magnitude_sum.is_finite() {
        return None;
    }
    let error_bound = (4.0 * f64::EPSILON * magnitude_sum).next_up();
    if determinant > error_bound {
        Some(1)
    } else if determinant < -error_bound {
        Some(-1)
    } else {
        None
    }
}

fn scalar_sign(value: f64) -> i8 {
    if value > 0.0 {
        1
    } else if value < 0.0 {
        -1
    } else {
        0
    }
}

fn exact_sum_of_products_sign<const N: usize>(terms: [(f64, f64, bool); N]) -> i8 {
    let products = terms.map(|(left, right, negate)| dyadic_product(left, right, negate));
    let Some(min_exponent) = products.iter().flatten().map(|term| term.exponent).min() else {
        return 0;
    };

    let mut sum = Integer::new();
    for product in products.into_iter().flatten() {
        let term = Integer::from(product.significand)
            << u32::try_from(product.exponent - min_exponent).unwrap();
        if product.negative {
            sum -= term;
        } else {
            sum += term;
        }
    }

    match sum.cmp0() {
        Ordering::Less => -1,
        Ordering::Equal => 0,
        Ordering::Greater => 1,
    }
}

pub fn f(value: impl Into<f64>) -> Float {
    Float::singleton(value)
}

pub fn cf(re: impl Into<f64>, im: impl Into<f64>) -> ComplexFloat {
    ComplexFloat::raw(f(re), f(im))
}

enum TrigFn {
    Sin,
    Cos,
}

fn trig_bounds(numerator: usize, denominator: usize, which: TrigFn) -> Float {
    assert!(denominator > 0, "zero denominator");
    let (arg_low, arg_high) = trig_argument_bounds(numerator, denominator);
    let (mut lower, mut upper) = match which {
        TrigFn::Sin => (
            crlibm::sinpi_rd(arg_low).min(crlibm::sinpi_rd(arg_high)),
            crlibm::sinpi_ru(arg_low).max(crlibm::sinpi_ru(arg_high)),
        ),
        TrigFn::Cos => (
            crlibm::cospi_rd(arg_low).min(crlibm::cospi_rd(arg_high)),
            crlibm::cospi_ru(arg_low).max(crlibm::cospi_ru(arg_high)),
        ),
    };

    match which {
        TrigFn::Sin => {
            if interval_contains_point(arg_low, arg_high, 0.5) {
                upper = 1.0;
            }
            if interval_contains_point(arg_low, arg_high, 1.5) {
                lower = -1.0;
            }
        }
        TrigFn::Cos => {
            if interval_contains_point(arg_low, arg_high, 0.0)
                || interval_contains_point(arg_low, arg_high, 2.0)
            {
                upper = 1.0;
            }
            if interval_contains_point(arg_low, arg_high, 1.0) {
                lower = -1.0;
            }
        }
    }

    Float::from_bounds(lower, upper)
}

fn trig_argument_bounds(numerator: usize, denominator: usize) -> (f64, f64) {
    let reduced = Float::int(numerator % denominator);
    let denominator = Float::int(denominator);
    (
        round_down(2.0 * reduced.lower / denominator.upper),
        round_up(2.0 * reduced.upper / denominator.lower),
    )
}

fn interval_contains_point(lower: f64, upper: f64, point: f64) -> bool {
    lower <= point && point <= upper
}

fn interval_mul(a_low: &f64, a_high: &f64, b_low: &f64, b_high: &f64) -> Float {
    let lower = [
        round_down(a_low * b_low),
        round_down(a_low * b_high),
        round_down(a_high * b_low),
        round_down(a_high * b_high),
    ].into_iter().reduce(f64::min).unwrap();

    let upper = [
        round_up(a_low * b_low),
        round_up(a_low * b_high),
        round_up(a_high * b_low),
        round_up(a_high * b_high),
    ].into_iter().reduce(f64::max).unwrap();

    Float::from_bounds(lower, upper)
}

fn interval_div(a_low: &f64, a_high: &f64, b_low: &f64, b_high: &f64) -> Float {
    let lower = [
        round_down(a_low / b_low),
        round_down(a_low / b_high),
        round_down(a_high / b_low),
        round_down(a_high / b_high),
    ].into_iter().reduce(f64::min).unwrap();

    let upper = [
        round_up(a_low / b_low),
        round_up(a_low / b_high),
        round_up(a_high / b_low),
        round_up(a_high / b_high),
    ].into_iter().reduce(f64::max).unwrap();

    Float::from_bounds(lower, upper)
}

fn round_down(value: f64) -> f64 {
    if value.is_nan() {
        value
    } else {
        value.next_down()
    }
}

fn round_up(value: f64) -> f64 {
    if value.is_nan() {
        value
    } else {
        value.next_up()
    }
}

#[cfg(test)]
mod bounded_tests {
    use super::{cf, f, Float};
    use crate::numbers::{Complex, Real};
    use rug::Rational;

    fn rational(value: f64) -> Rational {
        Rational::from_f64(value).unwrap()
    }

    fn reference_sign(value: Rational) -> i8 {
        match value.cmp0() {
            std::cmp::Ordering::Less => -1,
            std::cmp::Ordering::Equal => 0,
            std::cmp::Ordering::Greater => 1,
        }
    }

    fn next_finite(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let mut bits = *state;
        if bits & (0x7ff_u64 << 52) == 0x7ff_u64 << 52 {
            bits ^= 1_u64 << 52;
        }
        f64::from_bits(bits)
    }

    #[test]
    fn singleton_predicates_resolve_rounded_cancellation() {
        let one_up = 1.0_f64.next_up();
        let one_down = 1.0_f64.next_down();

        let cross_left = cf(one_up, 1.0);
        let cross_right = cf(1.0, one_down);
        assert_eq!(one_up * one_down - 1.0, 0.0);
        assert_eq!(super::ComplexFloat::cross_sign(&cross_left, &cross_right), 1);

        let dot_left = cf(one_up, 1.0);
        let dot_right = cf(one_down, -1.0);
        assert_eq!(one_up * one_down - 1.0, 0.0);
        assert_eq!(super::ComplexFloat::dot_sign(&dot_left, &dot_right), 1);

        let origin = cf(0.0, 0.0);
        let orientation_left = cf(1.0, 1.0);
        let orientation_right = cf(2.0, 2.0_f64.next_up());
        assert_eq!(
            super::ComplexFloat::orientation_sign(
                &origin,
                &orientation_left,
                &orientation_right,
            ),
            1,
        );

        assert_eq!(
            super::ComplexFloat::cross_sign(&cf(1.0, 1.0), &cf(2.0, 2.0)),
            0,
        );
        assert_eq!(
            super::ComplexFloat::dot_sign(&cf(1.0, 1.0), &cf(1.0, -1.0)),
            0,
        );
        assert_eq!(
            super::ComplexFloat::orientation_sign(
                &cf(-0.0, 0.0),
                &cf(1.0, 1.0),
                &cf(2.0, 2.0),
            ),
            0,
        );
    }

    #[test]
    fn orientation_filter_rejects_a_wrong_rounded_sign() {
        let offset = 2.0_f64.powi(-51);
        let origin = cf(offset, -offset);
        let left = cf(1.0, 2.0);
        let right = cf(2.0_f64.next_down(), 4.0_f64.next_up());
        let naive = (1.0 - offset) * (4.0_f64.next_up() + offset)
            - (2.0 + offset) * (2.0_f64.next_down() - offset);

        assert!(naive > 0.0);
        assert_eq!(
            super::ComplexFloat::orientation_sign(&origin, &left, &right),
            -1,
        );
    }

    #[test]
    fn singleton_predicates_handle_underflow_and_overflow() {
        let tiny = f64::from_bits(1);
        assert_eq!(
            super::ComplexFloat::cross_sign(&cf(tiny, 0.0), &cf(0.0, tiny)),
            1,
        );
        assert_eq!(
            super::ComplexFloat::dot_sign(&cf(tiny, tiny), &cf(tiny, tiny)),
            1,
        );
        assert_eq!(
            super::ComplexFloat::orientation_sign(
                &cf(0.0, 0.0),
                &cf(tiny, 0.0),
                &cf(0.0, tiny),
            ),
            1,
        );

        let maximum = f64::MAX;
        let below_maximum = maximum.next_down();
        assert_eq!(
            super::ComplexFloat::cross_sign(
                &cf(maximum, maximum),
                &cf(maximum, below_maximum),
            ),
            -1,
        );
        assert_eq!(
            super::ComplexFloat::dot_sign(
                &cf(maximum, maximum),
                &cf(maximum, -below_maximum),
            ),
            1,
        );
        assert_eq!(
            super::ComplexFloat::orientation_sign(
                &cf(maximum, 0.0),
                &cf(-maximum, 1.0),
                &cf(0.0, 2.0),
            ),
            -1,
        );
    }

    #[test]
    fn singleton_predicates_match_exact_rational_arithmetic() {
        let mut state = 0x243f_6a88_85a3_08d3;
        for _ in 0..10_000 {
            let values = std::array::from_fn::<_, 6, _>(|_| next_finite(&mut state));
            let [ox, oy, ax, ay, bx, by] = values;
            let origin = cf(ox, oy);
            let left = cf(ax, ay);
            let right = cf(bx, by);

            let expected_cross = reference_sign(
                rational(ax) * rational(by) - rational(ay) * rational(bx),
            );
            let expected_orientation = reference_sign(
                (rational(ax) - rational(ox)) * (rational(by) - rational(oy))
                    - (rational(ay) - rational(oy)) * (rational(bx) - rational(ox)),
            );
            let expected_dot = reference_sign(
                rational(ax) * rational(bx) + rational(ay) * rational(by),
            );

            assert_eq!(
                super::ComplexFloat::cross_sign(&left, &right),
                expected_cross,
                "cross mismatch for {values:?}",
            );
            assert_eq!(
                super::ComplexFloat::orientation_sign(&origin, &left, &right),
                expected_orientation,
                "orientation mismatch for {values:?}",
            );
            assert_eq!(
                super::ComplexFloat::dot_sign(&left, &right),
                expected_dot,
                "dot mismatch for {values:?}",
            );
        }
    }

    #[test]
    fn arithmetic_propagates_outward_bounds() {
        let a = Float::new(1.0, 1.25);
        let b = Float::new(-2.0, 3.0);
        let product = a * b;
        let (lower, upper) = product.bounds();
        assert!(lower <= -2.5);
        assert!(upper >= 3.75);
    }

    #[test]
    fn sign_is_unknown_when_interval_crosses_zero() {
        assert_eq!(Float::new(-1e-20, 1e-20).sign(), 0);
        assert_eq!(f(1.0).sign(), 1);
        assert_eq!(f(-1.0).sign(), -1);
    }

    #[test]
    fn trig_constructors_include_exact_extrema() {
        let cos_zero = Float::cos(0, 7);
        let (lower, upper) = cos_zero.bounds();
        assert!(lower <= 1.0 && upper >= 1.0);

        let sin_quarter = Float::sin(1, 4);
        let (lower, upper) = sin_quarter.bounds();
        assert!(lower <= 1.0 && upper >= 1.0);
    }

    #[test]
    fn integer_constructor_encloses_usize_max() {
        let (lower, upper) = Float::int(usize::MAX).bounds();
        let exact = Rational::from(usize::MAX);
        assert!(rational(lower) <= exact && exact <= rational(upper));
    }

    #[test]
    #[cfg(target_pointer_width = "64")]
    fn trig_constructors_enclose_large_integer_angles() {
        let numerator = 7_614_254_673_846_928_943usize;
        let denominator = 10_399_415_824_084_597_185usize;
        let argument = rug::Float::with_val(
            256,
            Rational::from((numerator, denominator)) * 2,
        );
        let expected_cos = argument.clone().cos_pi();
        let expected_sin = argument.sin_pi();
        for (enclosure, expected) in [
            (Float::cos(numerator, denominator), expected_cos),
            (Float::sin(numerator, denominator), expected_sin),
        ] {
            let (lower, upper) = enclosure.bounds();
            assert!(expected >= lower && expected <= upper, "{expected} outside [{lower}, {upper}]");
        }
    }

    #[test]
    fn complex_division_keeps_expected_value_in_bounds() {
        let quotient = cf(1.0, 1.0) / cf(1.0, -1.0);
        let (re_low, re_high) = quotient.re().bounds();
        let (im_low, im_high) = quotient.im().bounds();
        assert!(re_low <= 0.0 && re_high >= 0.0);
        assert!(im_low <= 1.0 && im_high >= 1.0);
    }
}
