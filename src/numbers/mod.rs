mod float;
mod exact;
use std::fmt::Debug;
use std::ops::{Add, Div, Mul, Sub, Neg};
pub use float::{Float, ComplexFloat, f, cf};
pub use exact::{ExactReal, ExactComplex};


pub trait TrigonometricHelper<R>: Clone + Debug {
	fn new(denominator: usize) -> Self;
	fn denominator(&self) -> usize;
	fn cos(&self, numerator: usize, denominator: usize) -> R;
	fn sin(&self, numerator: usize, denominator: usize) -> R;
}


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
		if self.sign() > 0 {
			self
		} else if self.sign() == 0 {
			Self::zero()
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

	fn abs2(&self) -> Self::Real {
		self.re().clone() * self.re().clone() + self.im().clone() * self.im().clone()
	}
	fn dir(angle_num: usize, angle_den: usize) -> Self {
		Self::raw(
			Self::Real::cos(angle_num, angle_den),
			Self::Real::sin(angle_num, angle_den)
		)
	}
	fn polar(r: Self::Real, dir: Self) -> Self {
		Self::from_real(r) * dir
	}
}

pub fn cross_product<C: Complex>(a: C, b: C) -> C::Real {
	a.re() * b.im() - a.im() * b.re()
}
