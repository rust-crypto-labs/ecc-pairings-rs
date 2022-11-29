use crate::{errors::ErrorKind, field::Field};

type WCoeffs<F> = (F, F, F, F, F, F);

// Generic elliptic curve
#[derive(Clone, Debug, PartialEq)]
pub struct EllipticCurve<F: Field> {
    weierstrass_coefficients: WCoeffs<F>,
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
    pub fn new_long_weierstrass(coeffs: WCoeffs<F>) -> Self {
        EllipticCurve {
            weierstrass_coefficients: coeffs,
        }
    }

    // Check that point is on the curve
    pub fn is_on_curve(self, p: &ECPoint<F>) -> bool {
        match p {
            ECPoint::PointAtInfinity => true,
            ECPoint::AffinePoint(x_p, y_p) => {
                let (a1, a2, a3, a4, _, a6) = self.get_a_invariants();
                // y² + a1 xy + a3 y = x³ + a2 x² + a4 x + a6
                y_p.square().add(&y_p.mul(&x_p.mul(a1))).add(&y_p.mul(a3))
                    == x_p
                        .zpow(3)
                        .add(&x_p.square().mul(a2))
                        .add(&x_p.mul(a4))
                        .add(a6)
            }
        }
    }

    // Random point
    pub fn random_point(self) -> ECPoint<F> {
        let (a1, a2, a3, a4, _, a6) = self.get_a_invariants();
        // Get a random x
        let rand_x = F::random_element();

        // y² + ( a1 x + a3 ) * y = x³ + a2 x² + a4 x + a6
        // b = a1 x + a3
        let b = &rand_x.mul(a1).add(a3);

        // c = - ( x³ + a2 x² + a4 x + a6 )
        let c = &rand_x
            .zpow(3)
            .add(&rand_x.square().mul(a2))
            .add(&rand_x.mul(a4))
            .add(a6)
            .neg();

        let delta = &b.square().add(&c.zmul(4).neg());

        let half = match c.one().zmul(2).invert() {
            Ok(x) => *x,
            Err(_) => todo!(),
        };

        // y = ( - b + sqrt( delta ) ) / 2
        let sq = match delta.sqrt() {
            Ok(x) => *x,
            Err(_) => todo!(),
        };
        let rand_y = half.mul(&b.neg().add(&sq));

        ECPoint::AffinePoint(rand_x, rand_y)
    }

    pub fn infinity_point() -> ECPoint<F> {
        ECPoint::PointAtInfinity
    }

    // Get long Weierstrass coeffs
    pub fn get_a_invariants(&self) -> &WCoeffs<F> {
        &self.weierstrass_coefficients
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

        let zero = x_r.zero();

        match (pt_p, pt_q) {
            // Case P = Q = 0
            (ECPoint::PointAtInfinity, ECPoint::PointAtInfinity) => Ok(zero),
            (ECPoint::PointAtInfinity, ECPoint::AffinePoint(x_q, _)) => {
                // Case P = 0
                // xR - xQ
                let x_q_neg = x_q.neg();
                Ok(x_r.add(&x_q_neg))
            }
            (ECPoint::AffinePoint(x_p, _), ECPoint::PointAtInfinity) => {
                // Case Q = 0
                // xR - xP
                let x_p_neg = x_p.neg();
                Ok(x_r.add(&x_p_neg))
            }
            (ECPoint::AffinePoint(x_p, y_p), ECPoint::AffinePoint(x_q, y_q)) => {
                let x_p_neg = x_p.neg();
                let y_p_neg = y_p.neg();

                if (x_p != x_q) || (y_p != y_q) {
                    // Case P != Q
                    if x_p == x_q {
                        // Case xP = xQ
                        // xR - xP
                        Ok(x_r.add(&x_p_neg))
                    } else {
                        // Case xP != xQ
                        let num = x_p_neg.add(y_q);
                        let denom = x_p_neg.add(x_q);
                        let slope = match num.div(&denom) {
                            Ok(x) => *x,
                            Err(_) => todo!(),
                        };

                        let xdiff = (x_r.add(&x_p_neg)).mul(&slope).neg();
                        let ydiff = y_r.add(&y_p_neg);
                        Ok(xdiff.add(&ydiff))
                    }
                } else {
                    // Case P = Q
                    let (a1, a2, a3, a4, _, _) = self.get_a_invariants();

                    // 3x² + 2x a2 - y a1 + a4
                    let num = x_p
                        .square()
                        .zmul(3)
                        .add(&y_p.mul(a1).neg())
                        .add(a4)
                        .add(&x_p.mul(a2).zmul(2));

                    // 2y + x a1 + a3
                    let denom = y_p.zmul(2).add(a3).add(&x_p.mul(a1));

                    if denom == denom.zero() {
                        // xR - xP
                        Ok(x_r.add(&x_p_neg))
                    } else {
                        let slope = match num.div(&denom) {
                            Ok(x) => *x,
                            Err(_) => todo!(),
                        };

                        let xdiff = (x_r.add(&x_p_neg)).mul(&slope).neg();
                        let ydiff = y_r.add(&y_p_neg);

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

        let (a1, a2, a3, a4, _, a6) = self.get_a_invariants();

        if x_p == x_q && (y_p.add(y_q).add(&a1.mul(x_q)).add(a3)).is_zero() {
            EllipticCurve::infinity_point()
        } else {
            let denom;
            let lambda_num;
            let nu_num;

            if x_p == x_q {
                denom = a3.add(&y_p.zmul(2)).add(&a1.mul(x_p));

                lambda_num = a4
                    .add(&x_p.square().zmul(3))
                    .add(&a2.mul(x_p).zmul(2))
                    .add(&a1.mul(y_p).neg());

                nu_num = (x_p.square().mul(x_p).neg())
                    .add(&a4.mul(x_p))
                    .add(&a6.zmul(2))
                    .add(&a3.mul(y_p).neg())
            } else {
                denom = x_q.add(&x_p.neg());
                lambda_num = y_q.add(&y_p.neg());
                nu_num = y_p.mul(x_q).add(&y_q.mul(x_p).neg());
            }

            let lambda = match lambda_num.div(&denom) {
                Ok(x) => *x,
                Err(_) => todo!(),
            };

            let nu = match nu_num.div(&denom) {
                Ok(x) => *x,
                Err(_) => todo!(),
            };

            let x = a2
                .neg()
                .add(&x_p.add(x_q).neg())
                .add(&lambda.square())
                .add(&a1.mul(&lambda));
            let y = a3.add(&nu).add(&x.mul(&lambda.add(a1))).neg();

            ECPoint::AffinePoint(x, y)
        }
    }

    // Doubles P
    pub fn double(&self, pt_p: &ECPoint<F>) -> ECPoint<F> {
        let (x_p, y_p) = match pt_p {
            ECPoint::PointAtInfinity => return pt_p.clone(),
            ECPoint::AffinePoint(x, y) => (x, y),
        };

        let (a1, a2, a3, a4, _, _) = self.get_a_invariants();

        let coeff = match y_p.zmul(2).add(&x_p.mul(a1)).add(a3).invert() {
            Ok(x) => *x,
            Err(_) => todo!(),
        };

        let lambda = &a1
            .zmul(3)
            .add(&x_p.mul(a2).zmul(2))
            .add(&y_p.mul(a1).neg())
            .add(a4)
            .mul(&coeff);

        let res_x = lambda
            .square()
            .add(&a1.mul(lambda))
            .add(&a2.neg())
            .add(&x_p.zmul(2).neg());

        let res_y = res_x
            .mul(a1)
            .neg()
            .add(&a3.neg())
            .add(&res_x.mul(lambda))
            .add(&x_p.mul(lambda))
            .add(&y_p.neg());

        ECPoint::AffinePoint(res_x, res_y)
    }

    // Returns the inverse of P
    // /!\ Can not invert zero point
    pub fn invert(&self, pt_p: &ECPoint<F>) -> Result<ECPoint<F>, ErrorKind> {
        let (x, y) = match pt_p {
            ECPoint::PointAtInfinity => return Err(ErrorKind::InvalidInput("P must not be zero")),
            ECPoint::AffinePoint(x, y) => (x, y),
        };
        let (a1, _, a3, _, _, _) = self.get_a_invariants();
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
}
