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

    // Power
    fn pow(self, y: &BigInt) -> Self;

    // Division
    fn div(self, y: &Self) -> Self;

    // Squaring
    fn square(self) -> Self;

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
}
