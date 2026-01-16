use std::ops::{Add, Div, Mul, Sub};
use ordered_float::OrderedFloat;

pub trait Real: Clone + PartialEq + Ord
+ Add<Output = Self> + Mul<Output = Self> + Sub<Output = Self> + Div<Output = Self> {
	fn zero() -> Self;
}

pub trait Complex: Clone + PartialEq
+ Add<Output = Self> + Mul<Output = Self> + Sub<Output = Self> + Div<Output = Self> {
	type Real: Real;
	fn re(&self) -> Self::Real;
	fn im(&self) -> Self::Real;
	fn raw(re: Self::Real, im: Self::Real) -> Self;
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

