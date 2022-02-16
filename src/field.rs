use rug::Integer;
/// Generic finite field operations
pub trait Field {
    /// Neutral element for addition
    fn zero() -> Self;

    /// Neutral element for multiplication
    fn one() -> Self;

    /// Addition
    fn add(self, y: &Self) -> Self;

    /// Multiplication
    fn mul(self, y: &Self) -> Self;

    /// Multiplication by an integer
    fn zmul(self, y: i64) -> Self;

    /// Power
    fn pow(self, y: &Integer) -> Self;

    /// Int power
    fn zpow(self, y: i64) -> Self;

    /// Division
    fn div(self, y: &Self) -> Self;

    /// Squaring
    fn square(self) -> Self;

    /// Square root
    fn sqrt(self) -> Result<Self, &'static str>;

    /// Multiplicative inverse
    fn invert(self) -> Self;

    /// Additive inverse
    fn neg(self) -> Self;

    /// Degree of the extension
    fn degree() -> usize;

    /// Field order
    fn order(self) -> u32;

    /// Base field order
    fn base_order() -> Integer;

    /// Random field point
    fn random_element() -> Self;
}

pub enum Scalar<const P:u32, const N:usize> {
    PFScalar(PrimeField<P>),
    FFScalar(FiniteField<P,N>),
}

#[derive(Clone, Debug, PartialEq)]
pub struct PrimeField<const P:u32> {
    pub value: Integer,
}

#[derive(Clone, Debug, PartialEq)]
pub struct FiniteField<const P: u32, const N: usize> {
    pub coords: [PrimeField<P>; N],
    pub polynomial: [PrimeField<P>; N],
}

impl<const P:u32> Field for PrimeField<P> {
    // Neutral element for addition
    fn zero(self) -> Self {
        let z = Integer::new();
        return PrimeField {value: z }
    }

    // Neutral element for multiplication
    fn one(self) -> Self {
        let one = Integer::from(1);
        return PrimeField {value: one};
    }

    // Addition
    fn add(self, y: &Self) -> Self {
        let res = self.value + y.value;
        PrimeField::<P> {
            value: res.modulo(&P),
        }
    }

    // Multiplication
    fn mul(&self, y: &Self) -> Self {
        let res = self.value * y.value;
        PrimeField::<P> {value: res.modulo(&P) }
    }

    // Multiplication by an integer
    fn zmul(self, y: &i64) -> Self {
        let res = self.value.mul(y.to_owned());
        PrimeField::<P> {value: res.modulo(&P) }
    }

    // Power
    fn pow(&self, y: &Integer) -> &Self {
        if y.is_zero() {
            return &self.one();
        } else if y.eq(&Integer::from(1)) {
            return &self;
        }

        if y.is_odd() {
            let n = (y.sub(&Integer::from(1))).div(2);
            return &self.pow(&n).mul(&self.pow(&n.add(Integer::from(1))));
        } else {
            let n = y.div(2);
            return &self.pow(&n).mul(&self.pow(&n));
        }
    }

    fn zpow(&self, y: i64) -> &Self {
        if y < 0 {
            return &self.invert().zpow(-y);
        }
        if y == 0 {
            return &self.one();
        } else if y==1 {
            return &self;
        } else if y%2 == 1 {
            let n = (y - 1) / 2;
            return &self.zpow(n).mul(&self.zpow(n+1));
        } else {
            let n = y / 2;
            return &self.zpow(n).mul(&self.zpow(n));
        }
    }

    // Division
    fn div(self, y: &Self) -> Self {
        return self.mul(&y.invert());
    }

    // Squaring
    fn square(self) -> Self {
        self.clone().mul(&self)
    }

    // Square root
    fn sqrt(self) -> Result<Self, &'static str> {
        if self.clone().zpow((i64::from(self.clone().order()) - 1)/2).to_owned() != self.clone().one() {
            return Err("Non quadratic residue.");
        }

        let one = self.clone().one();
        let q = self.clone().order() - 1;
        let s = 0;

        // Find Q, S such that p - 1 = Q * 2^S
        while q % 2 == 0 {
            s += 1;
            q /= 2;
        }

        let z = &one;
        while z.zpow((i64::from(self.clone().order()) - 1) / 2 ) == &one {
            z = &PrimeField::random_element();
        }
        let m = &s;
        let c  =z.zpow(i64::from(q));
        let t = self.clone().zpow(i64::from(q));
        let r = self.clone().zpow((i64::from(q) + 1)/2);
        let i = 0;
        while t.to_owned() != self.clone().one() && t.to_owned() != self.clone().zero() {
            let b = c.zpow(2^(m - i - 1));
            m = &i;
            c = &b.square();
            t = &t.mul(&c);
            r = &r.mul(&b);
        }

        if t.to_owned() == self.clone().zero() {
            return Ok(self.zero());
        }
        return Ok(r.zero());
    }

    // Multiplicative inverse
    fn invert(&self) -> Self {
        let ord = i64::from(&self.order() - 2);
        self.zpow(ord).to_owned()
    }

    // Additive inverse
    fn neg(self) -> Self {
        PrimeField { value: self.value.neg() }
    }

    // Degree of the extension
    fn degree() -> usize {
        return 1;
    }

    // Field order
    fn order(self) -> u32 {
        P
    }

    // Base field order
    fn base_order() -> Integer {
        P
    }

    // Random field point
    fn random_element() -> Self {
        todo!()
    }
}

impl<const P: u32, const N: usize> Field for FiniteField<P,N> {
    // Neutral element for addition
    fn zero(self) -> Self {
        let zero = PrimeField { value: Integer::new() };
        return FiniteField::<P,N> { coords: [zero; N], polynomial: self.polynomial }
    }

    // Neutral element for multiplication
    fn one(self) -> Self {
        let one = PrimeField { value: Integer::from(1) };
        let zero = PrimeField{value: Integer::new()};
        let res = [zero; N];
        res[0] = one;
        return FiniteField::<P,N> { coords: res, polynomial: self.polynomial }
    }

    fn add(self, y: &Self) -> Self {
        let x = self.clone().coords;
        let v = y.coords;
        for i in 0..N {
            x[i as usize] = x[i as usize].add(&v[i as usize]);
        }
        return FiniteField::<P,N> {coords: x, polynomial: self.polynomial};
    }

    fn mul(&self, y: &Self) -> Self {
        let zero = PrimeField { value: Integer::new() };
        // Initialize polynomials
        let A = self.clone().coords;
        let B = y.coords;
        let I = self.clone().polynomial;

        // Create a polynomial of degree 2N - 2
        // Must be a vec until const generic operations are allowed
        let Q = vec![zero; 2 * N - 1];

        // Create remainder polynomial
        let R = [zero; N];

        // Polynomial multiplication A * B
        for k in 0..(2 * N - 1) {
            for l in 0..(k+1) {
                if k - l < N && l < N {
                    Q[k] = Q[k].add(&A[k-l].mul(&B[l]));
                }
            }
        }

        // Polynomial euclidian remainder
        for l in N..(2*N-1) {
            let r = 2*N - 2 - l;
            for k in 0..N {
                Q[k+r-N] = Q[k+r-N].add(&Q[r].mul(&I[k]));
            }
        }

        for i  in 0..N {
            R[i] = Q[i];
        }
        FiniteField { coords:R, polynomial: self.polynomial }
    }

    fn zmul(self, y: &i64) -> Self {
        let x = self.clone().coords;
        for i in 0..N {
            x[i as usize] = x[i as usize].zmul(y)
        }
        return FiniteField::<P,N> {coords: x, polynomial: self.polynomial};
    }

    fn pow(&self, y: &Integer) -> &Self {
        if y.is_zero() {
            return &self.one();
        } else if y.eq(&Integer::from(1)) {
            return &self;
        }

        if y.is_odd() {
            let n = (y.sub(&Integer::from(1))).div(2);
            return &self.pow(&n).mul(&self.pow(&n.add(Integer::from(1))));
        } else {
            let n = y.div(2);
            return &self.pow(&n).mul(&self.pow(&n));
        }
    }

    fn zpow(&self, y: i64) -> &Self {
        if y == 0 {
            return &self.one();
        } else if y==1 {
            return &self;
        } else if y%2 == 1 {
            let n = (y - 1) / 2;
            return &self.zpow(n).mul(&self.zpow(n+1));
        } else {
            let n = y / 2;
            return &self.zpow(n).mul(&self.zpow(n));
        }
    }

    fn div(self, y: &Self) -> Self {
        return self.mul(&y.invert());
    }

    fn square(self) -> Self {
        self.clone().mul(&self)
    }

    fn sqrt(self) -> Result<Self, &'static str> {
        if self.clone().zpow((i64::from(self.clone().order()) - 1)/2).to_owned() != self.clone().one() {
            return Err("Non quadratic residue.");
        }

        let one = self.clone().one();
        let q = self.clone().order() - 1;
        let s = 0;

        // Find Q, S such that p - 1 = Q * 2^S
        while q % 2 == 0 {
            s += 1;
            q /= 2;
        }

        let z = &one;
        while z.zpow((i64::from(self.clone().order()) - 1) / 2 ) == &one {
            z = &FiniteField::random_element();
        }
        let m = &s;
        let c  =z.zpow(i64::from(q));
        let t = self.clone().zpow(i64::from(q));
        let r = self.clone().zpow((i64::from(q) + 1)/2);
        let i = 0;
        while t.to_owned() != self.clone().one() && t.to_owned() != self.clone().zero() {
            let b = c.zpow(2^(m - i - 1));
            m = &i;
            c = &b.square();
            t = &t.mul(&c);
            r = &r.mul(&b);
            i += 1;
        }

        if t.to_owned() == self.clone().zero() {
            return Ok(self.zero());
        }
        return Ok(r.zero());
    }

    fn invert(&self) -> Self {
        let ord = i64::from(&self.order() - 2);
        self.zpow(ord).to_owned()
    }

    fn neg(self) -> Self {
        let coords = self.clone().coords;
        for i in 0..N as usize {
            coords[i] = coords[i].neg();
        }

        return FiniteField {coords: coords, polynomial: self.polynomial};
    }


    fn degree() -> usize {
        N
    }

    fn order(self) -> u32 {
        todo!()
    }

    fn base_order() -> Integer {
        Integer.from(P)
    }

    fn random_element() -> Self {
        todo!()
    }
}