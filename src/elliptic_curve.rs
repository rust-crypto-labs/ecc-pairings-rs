use crate::field::Field;

// Generic elliptic curve
#[derive(Clone, Debug, PartialEq)]
pub struct EllipticCurve<F: Field> {
    weierstrass_coefficients: [F; 6],
}

// Rational point on an elliptic curve (affine coords)
#[derive(Clone, Debug, PartialEq)]
pub struct ECPoint<F: Field + Clone> {
    pub curve: EllipticCurve<F>,
    pub x: F,
    pub y: F,
}

// Elliptic curve data structure
impl<F: Field + Clone + PartialEq> EllipticCurve<F> {
    // New curve, long Weierstrass form
    // y² + a1 xy + a3 y = x³ + a2 x² + a4 x + a6
    pub fn new_long_weierstrass(coeffs: [F; 6]) -> Self {
        return EllipticCurve {
            weierstrass_coefficients: coeffs,
        };
    }

    // Random point
    pub fn random_point(self) -> ECPoint<F> {
        unimplemented!()
    }

    pub fn infinity_point(self) -> ECPoint<F> {
        unimplemented!()
    }

    // Get long Weierstrass coeffs
    pub fn get_a_invariants(&self) -> [F; 6] {
        return self.weierstrass_coefficients.clone();
    }
}

// Point operations
impl<F: Field + Clone + PartialEq> ECPoint<F> {
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
        return self.curve == other.curve && self.x == other.x && self.y == other.y;
    }

    // Returns the evaluation of the line PQ at R, where P is self
    // /!\ R cannot be the zero point
    pub fn line(&self, pt_q: &Self, pt_r: &Self) -> F {
        let x_p = self.x.clone();
        let y_p = self.y.clone();

        let x_q = pt_q.x.clone();
        let y_q = pt_q.y.clone();

        let x_r = pt_r.x.clone();
        let y_r = pt_r.y.clone();

        if self.is_zero() || pt_q.is_zero() {
            if self.is_equal(&pt_q) {
                // Case P = Q = 0
                // 1
                return F::one();
            } else if self.is_zero() {
                // Case P = 0
                // xR - xQ
                return x_r.clone().add(&x_q.clone().neg());
            } else {
                // Case Q = 0
                // xR - xP
                return x_r.clone().add(&x_p.clone().neg());
            }
        } else if !self.is_equal(&pt_q) {
            // Case P != Q
            if x_p == x_q {
                // Case xP = xQ
                // xR - xP
                return x_r.clone().add(&x_p.neg());
            } else {
                // Case xP != xQ
                let num = y_q.clone().add(&x_p.clone().neg());
                let denom = x_q.clone().add(&x_p.clone().neg());
                let slope = num.div(&denom);

                let xdiff = (x_r.clone().add(&x_p.clone().neg())).mul(&slope).neg();
                let ydiff = y_r.clone().add(&y_p.clone().neg());
                let res = xdiff.add(&ydiff);

                return res;
            }
        } else {
            // Case P = Q

            let a = self.curve.get_a_invariants();
            let (a1, a2, a3, a4) = (&a[0], &a[1], &a[2], &a[3]);

            // 3x² + 2x a2 - y a1 + a4
            let num = x_p
                .clone()
                .square()
                .zmul(3)
                .add(&y_p.clone().mul(&a1).neg())
                .add(&a4)
                .add(&x_p.clone().mul(&a2).zmul(2));

            // 2y + x a1 + a3
            let denom = y_p.clone().zmul(2).add(&a3).add(&x_p.clone().mul(&a1));

            if denom == F::zero() {
                // xR - xP
                return x_r.clone().add(&x_p.neg());
            } else {
                let slope = num.div(&denom);

                let xdiff = (x_r.clone().add(&x_p.clone().neg())).mul(&slope).neg();
                let ydiff = y_r.clone().add(&y_p.clone().neg());

                let res = ydiff.add(&xdiff);

                return res;
            }
        }
    }

    // Returns the addition of self with Q
    pub fn add(&self, pt_q: &Self) -> Self {
        let x1 = &self.x;
        let y1 = &self.y;
        let x2 = &pt_q.x;
        let y2 = &pt_q.y;
        let a = &self.curve.get_a_invariants();
        let (a1, a2, a3, a4, a6) = (&a[0], &a[1], &a[2], &a[3], &a[5]);

        } else {
            let mut lambda;
            let mut nu;
            if x1 == x2 {
                lambda = (3 * x1^2 + 2 * &a[1] * x1 + &a[3] - &a[0] * y1)/(2 * y1 + &a[0] * x1 + &a[2]);
                nu = (- x1^3 + &a[3] * x1 + 2 * &a[5] - &a[2] * y1)/(2 * y1 + &a[0] * x1 + &a[2]);
            } else {
                lambda = (y2 - y1)/(x2 - x1);
                nu = (y1 * x2 - y2 * x1)/(x2 - x1);
            }
            let x3 = lambda^2 + &a[0] * lambda - &a[1] - x1 - x2;
            return ECPoint {
                curve: *self.curve,
                x: x3,
                y: -(lambda + &a[0]) * x3 - nu - &a[2]
            }
        }
    }

    // Doubles self
    pub fn double(&self) -> Self {
        unimplemented!()
    }

    // Returns the inverse of self
    pub fn invert(&self) -> Self {
        let a_invariants = &self.curve.get_a_invariants();
        return ECPoint{
            curve: *self.curve,
            x: &self.x,
            y: -&self.y - &a_invariants[0] * &self.x - &a_invariants[2]
        }
    }
}
