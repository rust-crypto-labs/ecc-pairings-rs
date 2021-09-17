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
        let a = self.get_a_invariants();
        let (a1, a2, a3, a4, a6) = (&a[0], &a[1], &a[2], &a[3], &a[5]);
        // Get a random x
        let rand_x = F::random_element();

        // y² + ( a1 x + a3 ) * y = x³ + a2 x² + a4 x + a6
        // b = a1 x + a3
        let b = &rand_x.clone().mul(&a1.clone()).add(&a3.clone());

        // c = - ( x³ + a2 x² + a4 x + a6 )
        let c = &rand_x
            .clone()
            .zpow(3)
            .add(&rand_x.clone().square().mul(&a2.clone()))
            .add(&rand_x.clone().mul(&a4.clone()))
            .add(&a6.clone())
            .neg();
        let half = F::one().zmul(2).invert();
        let delta = &b.clone().square().add(&c.clone().zmul(4).neg());

        // y = ( - b + sqrt( delta ) ) / 2
        let rand_y = half.mul(&b.clone().neg().add(&delta.clone().sqrt()));

        return ECPoint {
            curve: self.clone(),
            x: rand_x.clone(),
            y: rand_y.clone(),
        };
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
        return ECPoint {
            curve: curve.clone(),
            x: x.clone(),
            y: y.clone(),
        };
    }

    // Returns true if and only if this is the point at infinity
    pub fn is_zero(&self) -> bool {
        return self.clone().is_equal(&self.curve.clone().infinity_point());
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
        if self.clone().is_equal(&pt_q.clone()) {
            return self.double();
        }
        let x1 = &self.x;
        let y1 = &self.y;
        let x2 = &pt_q.x;
        let y2 = &pt_q.y;
        let a = &self.curve.get_a_invariants();
        let (a1, a2, a3, a4, a6) = (&a[0], &a[1], &a[2], &a[3], &a[5]);

        if x1 == x2 && y1.clone().add(&y2).add(&a1.clone().mul(&x2)).add(&a3) == F::zero() {
            return self.curve.clone().infinity_point();
        } else {
            let lambda;
            let nu;
            if x1 == x2 {
                lambda = (a4
                    .clone()
                    .add(&x1.clone().square().zmul(3))
                    .add(&a2.clone().mul(&x1).zmul(2))
                    .add(&a1.clone().mul(&y1).neg()))
                .div(
                    &a3.clone()
                        .add(&y1.clone().zmul(2))
                        .add(&a1.clone().mul(&x1)),
                );
                nu = (x1.clone().square().mul(&x1).neg())
                    .add(&a4.clone().mul(&x1))
                    .add(&a6.clone().zmul(2))
                    .add(&a3.clone().mul(&y1).neg())
                    .div(
                        &a3.clone()
                            .add(&y1.clone().zmul(2))
                            .add(&a1.clone().mul(&x1)),
                    );
            } else {
                lambda = y2
                    .clone()
                    .add(&y1.clone().neg())
                    .div(&x2.clone().add(&x1.clone().neg()));
                nu = y1
                    .clone()
                    .mul(&x2)
                    .add(&y2.clone().mul(&x1).neg())
                    .div(&x2.clone().add(&x1.clone().neg()));
            }
            let x = a2
                .clone()
                .neg()
                .add(&x1.clone().add(&x2).neg())
                .add(&lambda.clone().square())
                .add(&a1.clone().mul(&lambda));
            let y = a3
                .clone()
                .add(&nu)
                .add(&x.clone().mul(&lambda.add(&a1)))
                .neg();

            return ECPoint {
                curve: self.curve.clone(),
                x,
                y,
            };
        }
    }

    // Doubles self
    pub fn double(&self) -> Self {
        let x_p = &self.x;
        let y_p = &self.y;
        let a = &self.curve.get_a_invariants();
        let (a1, a2, a3, a4) = (&a[0], &a[1], &a[2], &a[3]);

        let lambda = &a1
            .clone()
            .zmul(3)
            .add(&x_p.clone().mul(&a2.clone()).zmul(2))
            .add(&y_p.clone().mul(&a1.clone()).neg())
            .add(&a4.clone())
            .mul(
                &y_p.clone()
                    .zmul(2)
                    .add(&x_p.clone().mul(&a1.clone()))
                    .add(&a3.clone())
                    .invert(),
            );
        let res_x = &lambda
            .clone()
            .square()
            .add(&a1.clone().mul(&lambda.clone()))
            .add(&a2.clone().neg())
            .add(&x_p.clone().zmul(2).neg());
        return ECPoint {
            curve: self.curve.clone(),
            x: res_x.clone(),
            y: res_x
                .clone()
                .mul(&a1.clone())
                .neg()
                .add(&a3.clone().neg())
                .add(&res_x.clone().mul(&lambda.clone()))
                .add(&x_p.clone().mul(&lambda.clone()))
                .add(&y_p.clone().neg()),
        };
    }

    // Returns the inverse of self
    pub fn invert(&self) -> Self {
        let a = self.curve.get_a_invariants();
        let (a1, a3) = (&a[0], &a[2]);
        return ECPoint {
            curve: self.curve.clone(),
            x: self.x.clone(),
            y: a3.clone().add(&a1.clone().mul(&self.x)).add(&self.y).neg(),
        };
    }
}
