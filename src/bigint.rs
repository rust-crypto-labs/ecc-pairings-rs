use std::ops::Mul;
use std::ops::Add;

// Big number (to be replaced by an actual bignum library)
// Must support negative integers

#[derive(Clone, Debug, PartialEq)]
pub struct BigInt {}

// Example bignum operations
impl BigInt {
    // Returns the zero
    pub fn zero() -> Self {
        unimplemented!()
    }

    pub fn one() -> Self {
        unimplemented!()
    }

    // Returns the absolute value of self
    pub fn abs(&self) -> Self {
        unimplemented!()
    }

    // Returns true if the integer is zero
    pub fn is_zero(&self) -> bool {
        unimplemented!()
    }

    // Returns true if the integer is odd
    pub fn is_odd(&self) -> bool {
        unimplemented!()
    }

    // Returns true if the integer is positive
    pub fn is_positive(&self) -> bool {
        unimplemented!()
    }

    pub fn div(self, other: &Self) -> Self {
        unimplemented!()
    }

    pub fn pow(self, other: &Self) -> Self {
        unimplemented!()
    }

    pub fn sub(self, other: &Self) -> Self {
        unimplemented!()
    }

    // Returns the binary representation of self, starting
    // with the least-significant bits
    pub fn to_bits(self) -> Vec<bool> {
        unimplemented!()
    }

    pub fn modulo(self, y: &u32) -> Self {
        unimplemented!()
    }
}

impl Mul<BigInt> for BigInt {
    type Output = BigInt;
    fn mul(self, rhs: BigInt) -> Self::Output {
        unimplemented!()
    }
}

impl Mul<i64> for BigInt {
    type Output = BigInt;
    fn mul(self, rhs: i64) -> Self::Output {
        unimplemented!()
    }
}

impl Add<BigInt> for BigInt {
    type Output = BigInt;
    fn add(self, rhs: BigInt) -> Self::Output {
        unimplemented!()
    }
}