use crate::bigint::BigInt;
// Generic finite field operations
pub trait Field {
    // Neutral element for addition
    fn zero(self) -> Self;

    // Neutral element for multiplication
    fn one(self) -> Self;

    // Addition
    fn add(self, y: &Self) -> Self;

    // Multiplication
    fn mul(self, y: &Self) -> Self;

    // Multiplication by an integer
    fn zmul(self, y: &i64) -> Self;

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
    fn invert(&self) -> Self;

    // Additive inverse
    fn neg(self) -> Self;

    // Degree of the extension
    fn degree() -> u32;

    // Field order
    fn order(self) -> u32;

    // Base field order
    fn base_order() -> BigInt;

    // Random field point
    fn random_element() -> Self;
}

pub enum Scalar<const P:u32, const N:u32> {
    PFScalar(PrimeField<P>),
    FFScalar(FiniteField<P,N>),
}

#[derive(Clone, Debug, PartialEq)]
pub struct PrimeField<const P:u32> {
    pub value: BigInt,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteField<const P: u32, const N: u32> {
    pub coords: [PrimeField<P>; N],
    pub polynome: [PrimeField<P>; N],
}

impl<const P:u32> Field for PrimeField<P> {
    // Neutral element for addition
    fn zero(self) -> Self {
        return PrimeField::<P> {value: BigInt::zero() }
    }

    // Neutral element for multiplication
    fn one(self) -> Self {
        return PrimeField::<P> {value: BigInt::one()};
    }

    // Addition
    fn add(self, y: &Self) -> Self {
        let res = self.value + y.value;
        PrimeField::<P> {
            value: res.modulo(&P),
        }
    }

    // Multiplication
    fn mul(self, y: &Self) -> Self {
        let res = self.value * y.value;
        PrimeField::<P> {value: res.modulo(&P) }
    }

    // Multiplication by an integer
    fn zmul(self, y: &i64) -> Self {
        let res = self.value * y;
        PrimeField::<P> {value: res.modulo(&P) }
    }

    // Power
    fn pow(self, y: &BigInt) -> Self;

    // Int power
    fn zpow(self, y: i64) -> Self;

    // Division
    fn div(self, y: &Self) -> Self {
        return self.mul(&y.invert());
    }

    // Squaring
    fn square(self) -> Self;

    // Square root
    fn sqrt(self) -> Self;

    // Multiplicative inverse
    fn invert(&self) -> Self;

    // Additive inverse
    fn neg(self) -> Self;

    // Degree of the extension
    fn degree() -> u32 {
        return 1;
    }

    // Field order
    fn order(self) -> u32 {
        P
    }

    // Base field order
    fn base_order() -> BigInt {
        P
    }

    // Random field point
    fn random_element() -> Self;
}

impl<const P: u32, const N: u32> Field for FiniteField<P,N> {
    // Neutral element for addition
    fn zero(self) -> Self {
        let zero = PrimeField { value: BigInt::zero() };
        return FiniteField::<P,N> { coords: [zero; N], polynome: self.polynome }
    }

    // Neutral element for multiplication
    fn one(self) -> Self {
        let one = PrimeField { value: BigInt::one() };
        return FiniteField::<P,N> { coords: [one; N], polynome: self.polynome }
    }

    fn add(self, y: &Self) -> Self {
        let x = self.clone().coords;
        let v = y.coords;
        for i in 0..N {
            x[i as usize] = x[i as usize].add(&v[i as usize]);
        }
        return FiniteField::<P,N> {coords: x, polynome: self.polynome};
    }

    fn mul(self, y: &Self) -> Self {
        todo!()
    }

    fn zmul(self, y: &i64) -> Self {
        let x = self.clone().coords;
        for i in 0..N {
            x[i as usize] = x[i as usize].zmul(y)
        }
        return FiniteField::<P,N> {coords: x, polynome: self.polynome};
    }

    fn pow(self, y: &BigInt) -> Self {
        todo!()
    }

    fn zpow(self, y: i64) -> Self {
        todo!()
    }

    fn div(self, y: &Self) -> Self {
        todo!()
    }

    fn square(self) -> Self {
        todo!()
    }

    fn sqrt(self) -> Self {
        todo!()
    }

    fn invert(&self) -> Self {
        todo!()
    }

    fn neg(self) -> Self {
        todo!()
    }

    fn degree() -> u32 {
        N
    }

    fn order(self) -> u32 {
        todo!()
    }

    fn base_order() -> BigInt {
        P
    }

    fn random_element() -> Self {
        todo!()
    }
}