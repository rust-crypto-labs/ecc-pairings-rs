use crate::errors::ErrorKind;
use rug::{Integer, ops::DivRounding, rand::{RandState}};

/// Generic finite field operations
pub trait Field {
    /// Neutral element for addition
    fn zero(&self) -> Self;

    /// Neutral element for multiplication
    fn one(&self) -> Self;

    /// Check if value is zero
    fn is_zero(&self) -> bool;

    /// Addition
    fn add(&self, y: &Self) -> Self;

    /// Additive inverse
    fn neg(&self) -> Self;

    /// Multiplication
    fn mul(&self, y: &Self) -> Self;

    /// Multiplication by an integer
    fn zmul(&self, y: i64) -> Self;

    /// Power
    fn pow(&self, y: &Integer) -> Self;

    /// Int power
    fn zpow(&self, y: i64) -> Self;

    /// Squaring
    fn square(&self) -> Self;

    /// Square root
    fn sqrt(&self) -> Result<Box<Self>, ErrorKind>;

    fn is_square(&self) -> bool;

    /// Multiplicative inverse
    fn invert(&self) -> Result<Box<Self>, ErrorKind>;

    /// Division
    fn div(&self, y: &Self) -> Result<Box<Self>, ErrorKind>;

    /// Degree of the extension
    fn degree() -> usize;

    /// Field order
    fn order(&self) -> u32;

    /// Base field order
    fn base_order() -> Integer;

    /// Random field point
    fn random_element(&self) -> Self;
}

pub enum Scalar<const P: u32, const N: usize> {
    PFScalar(PrimeField<P>),
    FFScalar(FiniteField<P, N>),
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct PrimeField<const P: u32> {
    pub value: Integer,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteField<const P: u32, const N: usize> {
    pub coords: Vec<PrimeField<P>>,
    pub polynomial: Vec<PrimeField<P>>,
}

impl<const P: u32> Field for PrimeField<P> {
    // Neutral element for addition
    fn zero(&self) -> Self {
        let z = Integer::new();
        PrimeField { value: z }
    }

    // Neutral element for multiplication
    fn one(&self) -> Self {
        let one = Integer::from(1);
        PrimeField { value: one }
    }

    // Check if value is zero
    fn is_zero(&self) -> bool {
        self == &self.zero()
    }

    // Addition
    fn add(&self, y: &Self) -> Self {
        PrimeField::<P> {
            value: Integer::from(&self.value + &y.value) % P,
        }
    }

    // Multiplication
    fn mul(&self, y: &Self) -> Self {
        PrimeField::<P> {
            value: Integer::from(&self.value * &y.value) % P,
        }
    }

    // Multiplication by an integer
    fn zmul(&self, y: i64) -> Self {
        PrimeField::<P> {
            value: Integer::from(&self.value * y) % P,
        }
    }

    // Power
    fn pow(&self, y: &Integer) -> Self {
        let zero: u8 = 0;
        let one: u8 = 1;
        let two: u8 = 2;

        if y == &zero {
            return self.one();
        } else if y == &one {
            return self.clone();
        }

        let n = Integer::from(y / two);

        if y.is_odd() {
            self.pow(&n).mul(&self.pow(&(n + one)))
        } else {
            self.pow(&n).mul(&self.pow(&n))
        }
    }

    fn zpow(&self, y: i64) -> Self {
        if self.is_zero() {
            return self.zero();
        }

        if y < 0 {
            // Zero-cheking has been done prior, unwrap is valid
            self.invert().unwrap().zpow(-y)
        } else if y == 0 {
            self.one()
        } else if y == 1 {
            self.clone()
        } else if y % 2 == 1 {
            let n = (y - 1) / 2;
            self.zpow(n).mul(&self.zpow(n + 1))
        } else {
            let n = y / 2;
            self.zpow(n).mul(&self.zpow(n))
        }
    }

    // Division
    fn div(&self, y: &Self) -> Result<Box<Self>, ErrorKind> {
        y.invert().map(|x| Box::new(self.mul(&*x)))
    }

    // Squaring
    fn square(&self) -> Self {
        self.clone().mul(self)
    }

    // Square root
    fn sqrt(&self) -> Result<Box<Self>, ErrorKind> {
        let zero = self.zero();
        let one = self.one();
        let mut q = i64::from(self.order()) - 1;

        if self.zpow(q / 2) != one {
            return Err(ErrorKind::NonQuadraticResidue);
        }

        let mut s = 0;

        // Find Q, S such that p - 1 = Q * 2^S
        while q % 2 == 0 {
            s += 1;
            q /= 2;
        }

        let mut z = one.clone();
        while z.zpow(q / 2) == one {
            z = z.random_element();
        }
        let mut m = s;
        let mut c = z.zpow(q);
        let mut t = self.clone().zpow(q);
        let mut r = self.clone().zpow((q + 1) / 2);
        let mut i = 0;
        while t != one && t != zero {
            let b = c.zpow(2 ^ (m - i - 1));
            m = i;
            c = b.square();
            t = t.mul(&c);
            r = r.mul(&b);
            i += 1;
        }

        // Warning: zero ???
        if t.is_zero() {
            return Ok(Box::new(self.zero()));
        }
        Ok(Box::new(r.zero()))
    }

    fn is_square(&self) -> bool {
        // Legendre's symbol
        self.value.legendre(&Integer::from_f32(P as f32).unwrap()) == 1
    }

    // Multiplicative inverse
    fn invert(&self) -> Result<Box<Self>, ErrorKind> {
        if self.is_zero() {
            return Err(ErrorKind::NoInverse);
        }
        let ord = i64::from(&self.order() - 2);

        Ok(Box::new(self.zpow(ord)))
    }

    // Additive inverse
    fn neg(&self) -> Self {
        PrimeField {
            value: Integer::from(-&self.value),
        }
    }

    // Degree of the extension
    fn degree() -> usize {
        1
    }

    // Field order
    fn order(&self) -> u32 {
        P
    }

    // Base field order
    fn base_order() -> Integer {
        Integer::from(P)
    }

    // Random field point
    fn random_element(&self) -> Self {
        let rand = RandState::new();
        PrimeField {value: Integer::from(P).random_below(&mut rand )}
    }
}

impl<const P: u32, const N: usize> Field for FiniteField<P, N> {
    // Neutral element for addition
    fn zero(&self) -> Self {
        FiniteField::<P, N> {
            coords: Default::default(),
            polynomial: self.polynomial.clone(),
        }
    }

    // Neutral element for multiplication
    fn one(&self) -> Self {
        let one = PrimeField {
            value: Integer::from(1),
        };
        let mut res: Vec<PrimeField<P>> = vec![Default::default(); N];
        res[0] = one;

        FiniteField::<P, N> {
            coords: res,
            polynomial: self.polynomial.clone(),
        }
    }

    // Check if value is zero
    fn is_zero(&self) -> bool {
        self == &self.zero()
    }

    fn add(&self, y: &Self) -> Self {
        let mut x = self.coords.clone();
        let v = y.coords.clone();
        for i in 0..N {
            x[i] = x[i].add(&v[i]);
        }

        FiniteField::<P, N> {
            coords: x,
            polynomial: self.polynomial.clone(),
        }
    }

    fn mul(&self, y: &Self) -> Self {
        // Initialize polynomials
        let a = self.coords.clone();
        let b = y.coords.clone();
        let i = self.polynomial.clone();

        // Create a polynomial of degree 2N - 2
        // Must be a vec until const generic operations are allowed
        let mut q: Vec<PrimeField<P>> = vec![Default::default(); 2 * N - 1];

        // Create remainder polynomial
        let mut r: Vec<PrimeField<P>> = vec![Default::default(); N];

        // Polynomial multiplication A * B
        for k in 0..(2 * N - 1) {
            for l in 0..(k + 1) {
                if k - l < N && l < N {
                    q[k] = q[k].add(&a[k - l].mul(&b[l]));
                }
            }
        }

        // Polynomial euclidian remainder
        for l in N..(2 * N - 1) {
            let r = 2 * N - 2 - l;
            for k in 0..N {
                q[k + r - N] = q[k + r - N].add(&q[r].mul(&i[k]));
            }
        }

        r[..N].clone_from_slice(&q[..N]);

        FiniteField {
            coords: r,
            polynomial: self.polynomial.clone(),
        }
    }

    fn zmul(&self, y: i64) -> Self {
        let mut x = self.coords.clone();

        x.iter_mut().take(N).for_each(|i| *i = i.zmul(y));

        FiniteField::<P, N> {
            coords: x,
            polynomial: self.polynomial.clone(),
        }
    }

    fn pow(&self, y: &Integer) -> Self {
        let zero: u8 = 0;
        let one: u8 = 1;
        let two: u8 = 2;

        if y == &zero {
            return self.one();
        } else if y == &one {
            return self.clone();
        }

        let n = Integer::from(y / two);

        if y.is_odd() {
            self.pow(&n).mul(&self.pow(&(n + one)))
        } else {
            self.pow(&n).mul(&self.pow(&n))
        }
    }

    fn zpow(&self, y: i64) -> Self {
        if y == 0 {
            self.one()
        } else if y == 1 {
            self.clone()
        } else if y % 2 == 1 {
            let n = (y - 1) / 2;
            self.zpow(n).mul(&self.zpow(n + 1))
        } else {
            let n = y / 2;
            self.zpow(n).mul(&self.zpow(n))
        }
    }

    fn square(&self) -> Self {
        self.clone().mul(self)
    }

    fn sqrt(&self) -> Result<Box<Self>, ErrorKind> {
        let zero = self.one();
        let one = self.one();
        let mut q = i64::from(self.order()) - 1;

        if self.zpow(q / 2) != one {
            return Err(ErrorKind::NonQuadraticResidue);
        }

        let mut s = 0;

        // Find Q, S such that p - 1 = Q * 2^S
        while q % 2 == 0 {
            s += 1;
            q /= 2;
        }

        let mut z = one.clone();
        while z.zpow(q / 2) == one {
            z = z.mul(&z);
        }
        let mut m = s;
        let mut c = z.zpow(q);
        let mut t = self.clone().zpow(q);
        let mut r = self.clone().zpow((q + 1) / 2);
        let mut i = 0;
        while t != one && t != zero {
            let b = c.zpow(2 ^ (m - i - 1));
            m = i;
            c = b.square();
            t = t.mul(&c);
            r = r.mul(&b);
            i += 1;
        }

        // Warning: zero ???
        if t == zero {
            return Ok(Box::new(zero));
        }
        Ok(Box::new(r.zero()))
    }

    fn is_square(&self) -> bool {
        // Euler's criteria
        return self.pow(&Integer::from_f32((P - 1).div_euc(2) as f32).unwrap()).eq(&self.one());
    }

    fn invert(&self) -> Result<Box<Self>, ErrorKind> {
        if self.is_zero() {
            return Err(ErrorKind::NoInverse);
        }
        let ord = i64::from(&self.order() - 2);

        Ok(Box::new(self.zpow(ord)))
    }

    fn div(&self, y: &Self) -> Result<Box<Self>, ErrorKind> {
        y.invert().map(|x| Box::new(self.mul(&*x)))
    }

    fn neg(&self) -> Self {
        let mut coords = self.clone().coords;

        coords.iter_mut().take(N).for_each(|i| *i = i.neg());

        FiniteField {
            coords,
            polynomial: self.polynomial.clone(),
        }
    }

    fn degree() -> usize {
        N
    }

    fn order(&self) -> u32 {
        todo!()
    }

    fn base_order() -> Integer {
        Integer::from(P)
    }

    fn random_element(&self) -> Self {
        let z: PrimeField<P> = PrimeField {
            value: Default::default()
        };
        FiniteField {
            coords: 
            vec![z.random_element(); N],
            polynomial: self.polynomial,
        }
    }
}
