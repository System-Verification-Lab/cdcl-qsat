mod minkowski;
mod annular;
mod compute;
mod polygon;

pub use annular::{AngleInterval, AnnularSection, minkowski_product, minkowski_product_all};
use polygon::ConvexPolygon;
pub use compute::{check_constraint, poly_check_constraint, AngleFraction, ConstraintResult,
evaluate_constraint};
