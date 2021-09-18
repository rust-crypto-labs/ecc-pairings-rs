use crate::{errors::ErrorKind, field::Field};

// Generic elliptic curve
#[derive(Clone, Debug, PartialEq)]
pub struct EllipticCurve<F: Field> {
    weierstrass_coefficients: [F; 6],
}

// Rational point on an elliptic curve (affine coords)
#[derive(Clone, Debug, PartialEq)]
pub enum ECPoint<F: Field + Clone> {
    AffinePoint(F, F),
    PointAtInfinity,
}

// Elliptic curve data structure
impl<F: Field + Clone + PartialEq> EllipticCurve<F> {
    // New curve, long Weierstrass form
    // y² + a1 xy + a3 y = x³ + a2 x² + a4 x + a6
    pub fn new_long_weierstrass(coeffs: [F; 6]) -> Self {
        EllipticCurve {
            weierstrass_coefficients: coeffs,
        }
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

        ECPoint::AffinePoint(rand_x, rand_y)
    }

    pub fn infinity_point() -> ECPoint<F> {
        ECPoint::PointAtInfinity
    }

    // Get long Weierstrass coeffs
    pub fn get_a_invariants(&self) -> [F; 6] {
        self.weierstrass_coefficients.clone()
    }

    // Returns the evaluation of the line PQ at R, where P is self
    // /!\ R cannot be the zero point
    pub fn line(
        &self,
        pt_p: &ECPoint<F>,
        pt_q: &ECPoint<F>,
        pt_r: &ECPoint<F>,
    ) -> Result<F, ErrorKind> {
        let (x_r, y_r) = match pt_r {
            // Case P = Q = 0
            ECPoint::PointAtInfinity => {
                return Err(ErrorKind::InvalidInput("R cannot be the point at infinity"))
            }
            ECPoint::AffinePoint(x, y) => (x, y),
        };

        match (pt_p, pt_q) {
            // Case P = Q = 0
            (ECPoint::PointAtInfinity, ECPoint::PointAtInfinity) => Ok(F::one()),
            (ECPoint::PointAtInfinity, ECPoint::AffinePoint(x_q, _)) => {
                // Case P = 0
                // xR - xQ
                let x_q_neg = x_q.clone().neg();
                Ok(x_r.clone().add(&x_q_neg))
            }
            (ECPoint::AffinePoint(x_p, _), ECPoint::PointAtInfinity) => {
                // Case Q = 0
                // xR - xP
                let x_p_neg = x_p.clone().neg();
                Ok(x_r.clone().add(&x_p_neg))
            }
            (ECPoint::AffinePoint(x_p, y_p), ECPoint::AffinePoint(x_q, y_q)) => {
                let x_p_neg = x_p.clone().neg();
                let y_p_neg = y_p.clone().neg();

                if (x_p != x_q) || (y_p != y_q) {
                    // Case P != Q
                    if x_p == x_q {
                        // Case xP = xQ
                        // xR - xP
                        Ok(x_r.clone().add(&x_p_neg))
                    } else {
                        // Case xP != xQ
                        let num = x_p_neg.clone().add(y_q);
                        let denom = x_p_neg.clone().add(x_q);
                        let slope = num.div(&denom);

                        let xdiff = (x_r.clone().add(&x_p_neg)).mul(&slope).neg();
                        let ydiff = y_r.clone().add(&y_p_neg);
                        Ok(xdiff.add(&ydiff))
                    }
                } else {
                    // Case P = Q
                    let a = self.get_a_invariants();
                    let (a1, a2, a3, a4) = (&a[0], &a[1], &a[2], &a[3]);

                    // 3x² + 2x a2 - y a1 + a4
                    let num = x_p
                        .clone()
                        .square()
                        .zmul(3)
                        .add(&y_p.clone().mul(a1).neg())
                        .add(a4)
                        .add(&x_p.clone().mul(a2).zmul(2));

                    // 2y + x a1 + a3
                    let denom = y_p.clone().zmul(2).add(a3).add(&x_p.clone().mul(a1));

                    if denom == F::zero() {
                        // xR - xP
                        Ok(x_r.clone().add(&x_p_neg))
                    } else {
                        let slope = num.div(&denom);

                        let xdiff = (x_r.clone().add(&x_p_neg)).mul(&slope).neg();
                        let ydiff = y_r.clone().add(&y_p_neg);

                        Ok(ydiff.add(&xdiff))
                    }
                }
            }
        }
    }

    // Returns the addition of P with Q
    pub fn add(&self, pt_p: &ECPoint<F>, pt_q: &ECPoint<F>) -> ECPoint<F> {
        let (x_p, y_p) = match pt_p {
            ECPoint::PointAtInfinity => return pt_q.clone(),
            ECPoint::AffinePoint(x, y) => (x, y),
        };
        let (x_q, y_q) = match pt_q {
            ECPoint::PointAtInfinity => return pt_p.clone(),
            ECPoint::AffinePoint(x, y) => (x, y),
        };

        let a = self.get_a_invariants();
        let (a1, a2, a3, a4, a6) = (&a[0], &a[1], &a[2], &a[3], &a[5]);

        if x_p == x_q && y_p.clone().add(y_q).add(&a1.clone().mul(x_q)).add(a3) == F::zero() {
            EllipticCurve::infinity_point()
        } else {
            let lambda;
            let nu;
            if x_p == x_q {
                lambda = (a4
                    .clone()
                    .add(&x_p.clone().square().zmul(3))
                    .add(&a2.clone().mul(x_p).zmul(2))
                    .add(&a1.clone().mul(y_p).neg()))
                .div(
                    &a3.clone()
                        .add(&y_p.clone().zmul(2))
                        .add(&a1.clone().mul(x_p)),
                );
                nu = (x_p.clone().square().mul(x_p).neg())
                    .add(&a4.clone().mul(x_p))
                    .add(&a6.clone().zmul(2))
                    .add(&a3.clone().mul(y_p).neg())
                    .div(
                        &a3.clone()
                            .add(&y_p.clone().zmul(2))
                            .add(&a1.clone().mul(x_p)),
                    );
            } else {
                lambda = y_q
                    .clone()
                    .add(&y_p.clone().neg())
                    .div(&x_q.clone().add(&x_p.clone().neg()));
                nu = y_p
                    .clone()
                    .mul(x_q)
                    .add(&y_q.clone().mul(x_p).neg())
                    .div(&x_q.clone().add(&x_p.clone().neg()));
            }
            let x = a2
                .clone()
                .neg()
                .add(&x_p.clone().add(x_q).neg())
                .add(&lambda.clone().square())
                .add(&a1.clone().mul(&lambda));
            let y = a3
                .clone()
                .add(&nu)
                .add(&x.clone().mul(&lambda.add(a1)))
                .neg();

            ECPoint::AffinePoint(x, y)
        }
    }

    // Doubles P
    pub fn double(&self, pt_p: &ECPoint<F>) -> ECPoint<F> {
        let (x_p, y_p) = match pt_p {
            ECPoint::PointAtInfinity => return pt_p.clone(),
            ECPoint::AffinePoint(x, y) => (x, y),
        };

        let a = self.get_a_invariants();
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
        let res_x = lambda
            .clone()
            .square()
            .add(&a1.clone().mul(&lambda.clone()))
            .add(&a2.clone().neg())
            .add(&x_p.clone().zmul(2).neg());

        let res_y = res_x
            .clone()
            .mul(&a1.clone())
            .neg()
            .add(&a3.clone().neg())
            .add(&res_x.clone().mul(&lambda.clone()))
            .add(&x_p.clone().mul(lambda))
            .add(&y_p.clone().neg());

        ECPoint::AffinePoint(res_x, res_y)
    }

    // Returns the inverse of P
    // /!\ Can not invert zero point
    pub fn invert(&self, pt_p: &ECPoint<F>) -> Result<ECPoint<F>, ErrorKind> {
        let (x, y) = match pt_p {
            ECPoint::PointAtInfinity => return Err(ErrorKind::InvalidInput("P must not be zero")),
            ECPoint::AffinePoint(x, y) => (x, y),
        };
        let [a1, _, a3, _, _, _] = self.get_a_invariants();
        let new_y = a3.add(&a1.mul(x)).add(y).neg();

        Ok(ECPoint::AffinePoint(x.clone(), new_y))
    }
}

// Point on a curve
impl<F: Field + Clone + PartialEq> ECPoint<F> {
    // New point from affine coords
    pub fn new_affine(x: F, y: F) -> Self {
        ECPoint::AffinePoint(x, y)
    }

    // Returns true if and only if this is the point at infinity
    pub fn is_zero(&self) -> bool {
        ECPoint::PointAtInfinity.is_equal(self)
    }

    // Returns true if the two points are equal
    pub fn is_equal(&self, other: &Self) -> bool {
        self == other
    }
}
