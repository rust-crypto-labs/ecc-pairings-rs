// Big number (to be replaced by an actual bignum library)
// Must support negative integers
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
}
