use crate::field::Field;
// Generic elliptic curve
#[derive(Clone, Debug)]
pub struct EllipticCurve<F: Field> {
    weierstrass_coefficients: [F; 6],
}

// Rational point on an elliptic curve (affine coords)
#[derive(Clone, Debug)]
pub struct ECPoint<F: Field> {
    pub curve: EllipticCurve<F>,
    pub x: F,
    pub y: F,
}

// Elliptic curve data structure
impl<F: Field> EllipticCurve<F> {
    // New curve, long Weierstrass form
    pub fn new_long_weierstrass(coeffs: [F; 6]) -> Self {
        unimplemented!()
    }

    // Random point
    pub fn random_point(self) -> ECPoint<F> {
        unimplemented!()
    }

    // Get long Weierstrass coeffs
    pub fn get_long_weierstrass(&self) -> [F; 6] {
        unimplemented!()
    }
}

// Point operations
impl<F: Field> ECPoint<F> {
    // New point from affine coords
    pub fn new_affine(curve: &EllipticCurve<F>, x: F, y: F) -> Self {
        unimplemented!()
    }

    // Returns true if and only if this is the point at infinity
    pub fn is_zero(&self) -> bool {
        unimplemented!()
    }

    // Returns true if the two points are equal
    pub fn is_equal(&self, other: &Self) -> bool {
        unimplemented!()
    }

    // Returns the evaluation of the line PQ at R, where P is self
    pub fn line(&self, pt_q: &Self, pt_r: &Self) -> F {
        unimplemented!()
    }

    // Returns the addition of self with Q
    pub fn add(&self, pt_q: &Self) -> Self {
        unimplemented!()
    }

    // Doubles self
    pub fn double(&self) -> Self {
        unimplemented!()
    }

    // Returns the inverse of self
    pub fn invert(&self) -> Self {
        unimplemented!()
    }
}
