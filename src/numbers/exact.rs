use std::cmp::Ordering;
use std::f64::consts::PI;
use std::ops::{Add, Div, Mul, Neg, Sub};

use rug::float::Round;
use rug::ops::{AddAssignRound, DivAssignRound, MulAssignRound};
use rug::{Float, Integer, Rational};

use super::{Complex, Real};

#[derive(Clone, PartialEq, Eq, Debug)]
struct Cyclotomic<const K: usize> {
    coeffs: Vec<Rational>,
}

impl<const K: usize> Cyclotomic<K> {
    fn degree() -> usize {
        assert!(K >= 1, "K must be at least 1");
        1usize << (K - 1)
    }

    fn order() -> usize {
        1usize << K
    }

    fn zero() -> Self {
        Self {
            coeffs: vec![Rational::from(0); Self::degree()],
        }
    }

    fn from_int(value: i128) -> Self {
        let mut coeffs = vec![Rational::from(0); Self::degree()];
        coeffs[0] = Rational::from(value);
        Self { coeffs }
    }

    fn zeta_pow(exp: usize) -> Self {
        let n = Self::degree();
        let order = Self::order();
        let reduced = exp % order;
        let (sign, idx) = if reduced >= n {
            (-1, reduced - n)
        } else {
            (1, reduced)
        };
        let mut coeffs = vec![Rational::from(0); n];
        coeffs[idx] = Rational::from(sign);
        Self { coeffs }
    }

    fn add(&self, other: &Self) -> Self {
        let mut coeffs = Vec::with_capacity(self.coeffs.len());
        for (a, b) in self.coeffs.iter().zip(other.coeffs.iter()) {
            coeffs.push(a.clone() + b.clone());
        }
        Self { coeffs }
    }

    fn sub(&self, other: &Self) -> Self {
        let mut coeffs = Vec::with_capacity(self.coeffs.len());
        for (a, b) in self.coeffs.iter().zip(other.coeffs.iter()) {
            coeffs.push(a.clone() - b.clone());
        }
        Self { coeffs }
    }

    fn neg(&self) -> Self {
        self.scale(Rational::from(-1))
    }

    fn scale(&self, value: Rational) -> Self {
        let mut coeffs = Vec::with_capacity(self.coeffs.len());
        for coeff in self.coeffs.iter() {
            coeffs.push(coeff.clone() * value.clone());
        }
        Self { coeffs }
    }

    fn mul(&self, other: &Self) -> Self {
        let n = Self::degree();
        let mut coeffs = vec![Rational::from(0); n];
        for i in 0..n {
            for j in 0..n {
                let mut idx = i + j;
                let sign = if idx >= n {
                    idx -= n;
                    -1
                } else {
                    1
                };
                let term = self.coeffs[i].clone() * other.coeffs[j].clone();
                coeffs[idx] = if sign == 1 {
                    coeffs[idx].clone() + term
                } else {
                    coeffs[idx].clone() - term
                };
            }
        }
        Self { coeffs }
    }

    fn conjugate(&self) -> Self {
        let n = Self::degree();
        let mut coeffs = vec![Rational::from(0); n];
        coeffs[0] = self.coeffs[0].clone();
        for i in 1..n {
            let target = n - i;
            coeffs[target] = coeffs[target].clone() - self.coeffs[i].clone();
        }
        Self { coeffs }
    }

    fn i_unit() -> Self {
        assert!(K >= 2, "i requires K >= 2");
        Self::zeta_pow(1usize << (K - 2))
    }

    fn inverse(&self) -> Self {
        let n = Self::degree();
        let mut mat = vec![vec![Rational::from(0); 2 * n]; n];
        for col in 0..n {
            let basis = Self::zeta_pow(col);
            let product = self.mul(&basis);
            for row in 0..n {
                mat[row][col] = product.coeffs[row].clone();
            }
        }
        for row in 0..n {
            mat[row][n + row] = Rational::from(1);
        }

        for pivot in 0..n {
            let mut pivot_row = None;
            for row in pivot..n {
                if !mat[row][pivot].is_zero() {
                    pivot_row = Some(row);
                    break;
                }
            }
            let Some(row) = pivot_row else {
                panic!("non-invertible element");
            };
            if row != pivot {
                mat.swap(row, pivot);
            }
            let pivot_val = mat[pivot][pivot].clone();
            for col in pivot..2 * n {
                mat[pivot][col] = mat[pivot][col].clone() / pivot_val.clone();
            }
            for row in 0..n {
                if row == pivot {
                    continue;
                }
                let factor = mat[row][pivot].clone();
                if factor.is_zero() {
                    continue;
                }
                for col in pivot..2 * n {
                    mat[row][col] =
                        mat[row][col].clone() - factor.clone() * mat[pivot][col].clone();
                }
            }
        }

        let mut coeffs = Vec::with_capacity(n);
        for row in 0..n {
            coeffs.push(mat[row][n].clone());
        }
        Self { coeffs }
    }

    fn to_f64_pair(&self) -> (f64, f64) {
        let order = Self::order() as f64;
        let step = 2.0 * PI / order;
        let mut re = 0.0;
        let mut im = 0.0;
        for (i, coeff) in self.coeffs.iter().enumerate() {
            if coeff.is_zero() {
                continue;
            }
            let angle = step * i as f64;
            re += coeff.to_f64() * angle.cos();
            im += coeff.to_f64() * angle.sin();
        }
        (re, im)
    }

    fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|coeff| coeff.is_zero())
    }

    fn real_interval(&self, prec: u32) -> (Float, Float) {
        let mut lower = Float::with_val(prec, 0);
        let mut upper = Float::with_val(prec, 0);
        let denom = Integer::from(1) << (K - 1);

        for (i, coeff) in self.coeffs.iter().enumerate() {
            if coeff.is_zero() {
                continue;
            }

            let coeff_low = Float::with_val_round(prec, coeff, Round::Down).0;
            let coeff_high = Float::with_val_round(prec, coeff, Round::Up).0;
            let mut angle = Float::with_val(prec, Integer::from(i));
            angle.div_assign_round(&denom, Round::Nearest);

            let mut cos_low = angle.clone();
            cos_low.cos_pi_round(Round::Down);
            let mut cos_high = angle;
            cos_high.cos_pi_round(Round::Up);

            let (term_low, term_high) = interval_mul(&coeff_low, &coeff_high, &cos_low, &cos_high);
            lower.add_assign_round(term_low, Round::Down);
            upper.add_assign_round(term_high, Round::Up);
        }

        (lower, upper)
    }
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct ExactReal<const K: usize>(Cyclotomic<K>);

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct ExactComplex<const K: usize>(Cyclotomic<K>);

impl<const K: usize> Add for ExactReal<K> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        ExactReal(self.0.add(&other.0))
    }
}

impl<const K: usize> Sub for ExactReal<K> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        ExactReal(self.0.sub(&other.0))
    }
}

impl<const K: usize> Mul for ExactReal<K> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        ExactReal(self.0.mul(&other.0))
    }
}

impl<const K: usize> Div for ExactReal<K> {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        ExactReal(self.0.mul(&other.0.inverse()))
    }
}

impl<const K: usize> Neg for ExactReal<K> {
    type Output = Self;
    fn neg(self) -> Self {
        ExactReal(self.0.neg())
    }
}

impl<const K: usize> PartialOrd for ExactReal<K> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<const K: usize> Ord for ExactReal<K> {
    fn cmp(&self, other: &Self) -> Ordering {
        let diff = self.clone() - other.clone();
        if diff.0.is_zero() {
            return Ordering::Equal;
        }

        let mut prec = 128u32.max((K as u32) + 16);
        loop {
            let (low, high) = diff.0.real_interval(prec);
            let zero = Float::with_val(prec, 0);
            if high < zero {
                return Ordering::Less;
            }
            if low > zero {
                return Ordering::Greater;
            }
            prec = prec.saturating_mul(2);
        }
    }
}

impl<const K: usize> Real for ExactReal<K> {
    fn zero() -> Self {
        ExactReal(Cyclotomic::zero())
    }

    fn int(from: usize) -> Self {
        ExactReal(Cyclotomic::from_int(from as i128))
    }

    fn cos(numerator: usize, denominator: usize) -> Self {
        assert!(denominator > 0, "zero denominator");
        assert!(
            denominator.is_power_of_two(),
            "denominator must be a power of two"
        );
        let order = Cyclotomic::<K>::order();
        assert!(order % denominator == 0, "denominator must divide 2^K");
        let step = order / denominator;
        let exp = (numerator % denominator) * step;
        let zeta = Cyclotomic::<K>::zeta_pow(exp);
        let zeta_conj = Cyclotomic::<K>::zeta_pow((order - exp) % order);
        ExactReal(zeta.add(&zeta_conj).scale(Rational::from((1, 2))))
    }

    fn sin(numerator: usize, denominator: usize) -> Self {
        assert!(denominator > 0, "zero denominator");
        assert!(
            denominator.is_power_of_two(),
            "denominator must be a power of two"
        );
        let order = Cyclotomic::<K>::order();
        assert!(order % denominator == 0, "denominator must divide 2^K");
        let step = order / denominator;
        let exp = (numerator % denominator) * step;
        let zeta = Cyclotomic::<K>::zeta_pow(exp);
        let zeta_conj = Cyclotomic::<K>::zeta_pow((order - exp) % order);
        let diff = zeta.sub(&zeta_conj);
        let neg_i = Cyclotomic::<K>::i_unit().neg();
        ExactReal(diff.mul(&neg_i).scale(Rational::from((1, 2))))
    }

    fn to_float(&self) -> f64 {
        self.0.to_f64_pair().0
    }
}

impl<const K: usize> Add for ExactComplex<K> {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        ExactComplex(self.0.add(&other.0))
    }
}

impl<const K: usize> Sub for ExactComplex<K> {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        ExactComplex(self.0.sub(&other.0))
    }
}

impl<const K: usize> Mul for ExactComplex<K> {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        ExactComplex(self.0.mul(&other.0))
    }
}

impl<const K: usize> Div for ExactComplex<K> {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        ExactComplex(self.0.mul(&other.0.inverse()))
    }
}

impl<const K: usize> Neg for ExactComplex<K> {
    type Output = Self;
    fn neg(self) -> Self {
        ExactComplex(self.0.neg())
    }
}

impl<const K: usize> Complex for ExactComplex<K> {
    type Real = ExactReal<K>;
    fn re(&self) -> Self::Real {
        let conj = self.0.conjugate();
        let sum = self.0.add(&conj).scale(Rational::from((1, 2)));
        ExactReal(sum)
    }

    fn im(&self) -> Self::Real {
        let conj = self.0.conjugate();
        let diff = self.0.sub(&conj);
        let neg_i = Cyclotomic::<K>::i_unit().neg();
        ExactReal(diff.mul(&neg_i).scale(Rational::from((1, 2))))
    }

    fn raw(re: Self::Real, im: Self::Real) -> Self {
        let i_unit = Cyclotomic::<K>::i_unit();
        ExactComplex(re.0.add(&im.0.mul(&i_unit)))
    }
}

fn interval_mul(a_low: &Float, a_high: &Float, b_low: &Float, b_high: &Float) -> (Float, Float) {
    let products_low = [
        mul_round(a_low, b_low, Round::Down),
        mul_round(a_low, b_high, Round::Down),
        mul_round(a_high, b_low, Round::Down),
        mul_round(a_high, b_high, Round::Down),
    ];
    let mut min = products_low[0].clone();
    for value in products_low.into_iter().skip(1) {
        if value < min {
            min = value;
        }
    }

    let products_high = [
        mul_round(a_low, b_low, Round::Up),
        mul_round(a_low, b_high, Round::Up),
        mul_round(a_high, b_low, Round::Up),
        mul_round(a_high, b_high, Round::Up),
    ];
    let mut max = products_high[0].clone();
    for value in products_high.into_iter().skip(1) {
        if value > max {
            max = value;
        }
    }

    (min, max)
}

fn mul_round(a: &Float, b: &Float, round: Round) -> Float {
    let mut prod = a.clone();
    prod.mul_assign_round(b, round);
    prod
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_real_arithmetic_preserves_large_rational_coefficients() {
        type R = ExactReal<2>;
        let power = R::int(1 << 32);
        let denominator = power.clone() * power;
        let small = R::int(1) / denominator.clone();
        let sum = small.clone() + small.clone();
        assert_eq!(sum, R::int(2) / denominator.clone());
        assert!(small > R::zero());
        let squared = small.clone() * small.clone();
        assert!(squared > R::zero());
        assert_eq!(squared * denominator.clone() * denominator.clone(), R::int(1));
        assert_eq!(small * denominator, R::int(1));
    }
}
