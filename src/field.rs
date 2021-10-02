use crate::bigint::BigInt;
// Generic finite field operations
pub trait Field {
    // Neutral element for addition
    fn zero() -> Self;

    // Neutral element for multiplication
    fn one() -> Self;

    // Addition
    fn add(self, y: &Self) -> Self;

    // Multiplication
    fn mul(self, y: &Self) -> Self;

    // Multiplication by an integer
    fn zmul(self, y: i64) -> Self;

    // Power
    fn pow(self, y: &BigInt) -> Self;

    // Int power
    fn zpow(self, y: i64) -> Self;

    // Division
    fn div(self, y: &Self) -> Self;

    // Squaring
    fn square(self) -> Self;

    // Square root
    fn sqrt(self) -> Self;

    // Multiplicative inverse
    fn invert(self) -> Self;

    // Additive inverse
    fn neg(self) -> Self;

    // Degree of the extension
    fn degree() -> Self;

    // Field order
    fn order() -> Self;

    // Base field order
    fn base_order() -> BigInt;

    // Random field point
    fn random_element() -> Self;
}

pub enum Scalar<F:Field + Clone> {
    PFScalar(PrimeField),
    FFScalar(FiniteField),
}
pub struct PrimeField {
    pub order: i32,
    pub value: BigInt,
}

pub struct FiniteField {
    pub order: i32,
    pub degree: i32,
    pub coords: Vec<PrimeField>,
}

impl Field for PrimeField {
    // Neutral element for addition
    fn zero() -> Self {
        return PrimeField { order: Self::order(), value: BigInt::zero() }
    }

    // Neutral element for multiplication
    fn one() -> Self {
        return PrimeField {order: Self::order(), value: BigInt::one()};
    }

    // Addition
    fn add(self, y: &Self) -> Self {
        let res = self.value + y.value;
        PrimeField {
            order: self.order.clone(),
            value: res.modulo(&self.order),
        }
    }

    // Multiplication
    fn mul(self, y: &Self) -> Self {
        let res = self.value * y.value;
        PrimeField { order: self.order.clone(), value: res.modulo(&self.order) }
    }

    // Multiplication by an integer
    fn zmul(self, y: i64) -> Self;

    // Power
    fn pow(self, y: &BigInt) -> Self;

    // Int power
    fn zpow(self, y: i64) -> Self;

    // Division
    fn div(self, y: &Self) -> Self;

    // Squaring
    fn square(self) -> Self;

    // Square root
    fn sqrt(self) -> Self;

    // Multiplicative inverse
    fn invert(self) -> Self;

    // Additive inverse
    fn neg(self) -> Self;

    // Degree of the extension
    fn degree() -> Self {
        return 1;
    }

    // Field order
    fn order(self) -> Self {
        return self.order;
    }

    // Base field order
    fn base_order() -> BigInt;

    // Random field point
    fn random_element() -> Self;
}