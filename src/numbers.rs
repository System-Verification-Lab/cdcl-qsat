use std::ops::{Add, Div, Mul, Sub, Neg};
use ordered_float::OrderedFloat;
use std::f64::consts::PI;
use std::fmt::Debug;

pub trait Real: Clone + Debug + PartialEq + Ord
+ Add<Output = Self> + Mul<Output = Self> + Sub<Output = Self>
+ Div<Output = Self> + Neg<Output = Self> {
	fn zero() -> Self;
	fn int(from: usize) -> Self;
	fn sign(&self) -> i8 {
		if *self > Self::zero() {
			1
		} else if *self < Self::zero() {
			-1
		} else {
			0
		}
	}
	fn abs(self) -> Self {
		if self.sign() >= 0 {
			self
		} else {
			-self
		}
	}
	fn cos(numerator: usize, denominator: usize) -> Self;
	fn sin(numerator: usize, denominator: usize) -> Self;
	fn to_float(&self) -> f64;
}

pub trait Complex: Clone + Debug + PartialEq
+ Add<Output = Self> + Mul<Output = Self> + Sub<Output = Self>
+ Div<Output = Self> + Neg<Output = Self>
{
	type Real: Real;
	fn re(&self) -> Self::Real;
	fn im(&self) -> Self::Real;
	fn raw(re: Self::Real, im: Self::Real) -> Self;
	fn zero() -> Self {
		Self::raw(Self::Real::zero(), Self::Real::zero())
	}

	fn from_real(r: Self::Real) -> Self {
		Self::raw(r, Self::Real::zero())
	}
	fn to_float_pair(&self) -> (f64, f64) {
		(self.re().to_float(), self.im().to_float())
	}
}

pub fn cross_product<C: Complex>(a: C, b: C) -> C::Real {
	a.re() * b.im() - a.im() * b.re()
}


// Real float
pub type Float = OrderedFloat<f64>;
impl Real for Float {
	fn zero() -> Self {
		OrderedFloat(0.0)
	}
	fn int(from: usize) -> Self {
		OrderedFloat(from as f64)
	}
	fn cos(numerator: usize, denominator: usize) -> Self {
		OrderedFloat((2.0 * PI * numerator as f64 / denominator as f64).cos())
	}
	fn sin(numerator: usize, denominator: usize) -> Self {
		OrderedFloat((2.0 * PI * numerator as f64 / denominator as f64).sin())
	}
	fn to_float(&self) -> f64 {
		self.0
	}
}

// Commplex float
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
}

// Helpers
pub fn f(value: impl Into<f64>) -> Float {
	Float::from(value.into())
}

pub fn cf(re: impl Into<f64>, im: impl Into<f64>) -> ComplexFloat {
	ComplexFloat::raw(f(re), f(im))
}