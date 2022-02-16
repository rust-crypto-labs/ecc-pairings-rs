use rug::{Complete, Integer};

use std::ops::{Div, Sub};

use crate::{
    elliptic_curve::{ECPoint, EllipticCurve},
    errors::ErrorKind,
    field::Field,
};

trait IntegerExt {
    fn is_zero(&self) -> bool;
    fn is_positive(&self) -> bool;
    fn large_pow(&self, other: &Self) -> Self;
    fn to_bits(self) -> Vec<bool>;
}

impl IntegerExt for Integer {
    fn is_zero(&self) -> bool {
        *self == 0
    }

    fn is_positive(&self) -> bool {
        let zero: Self = 0.into();
        *self > zero
    }

    fn large_pow(&self, _other: &Self) -> Self {
        unimplemented!()
    }

    fn to_bits(self) -> Vec<bool> {
        unimplemented!()
    }
}

/// Miller's algorithm
/// Returns f_{n,P}(Q) where div(f_{n,P}) = n(P) - ([n]P) - (n-1)(0)
pub fn miller<F: Field + Clone + PartialEq>(
    curve: &EllipticCurve<F>,
    pt_p: &ECPoint<F>,
    pt_q: &ECPoint<F>,
    n: &Integer,
) -> Result<F, ErrorKind> {
    // Basic checks
    if pt_p == &ECPoint::PointAtInfinity {
        return Err(ErrorKind::InvalidInput("P must not be zero"));
    }
    if pt_q == &ECPoint::PointAtInfinity {
        return Err(ErrorKind::InvalidInput("Q must not be zero"));
    }
    if n.is_zero() {
        return Ok(F::one());
    }

    // Negative values of n are allowed, in which case
    // Q is evaluated instead at (v_{[n]P} f_{n,P)})^(-1)
    let sign = n.is_positive();
    //let n = n.abs();
    let nbits = n.abs_ref().complete().to_bits();

    let one = F::one();

    let mut t = one;
    let mut i: usize = nbits.len() - 1; // Will not underflow because n != 0
    let mut pt_v = pt_p.clone();

    if i != 0 {
        i -= 1;

        // Miller loop
        loop {
            let pt_s = curve.double(&pt_v);
            let ell = curve.line(&pt_v, &pt_v, pt_q)?;
            let vee = curve.line(&pt_s, &curve.invert(&pt_s)?, pt_q)?;
            t = t.square().mul(&ell.div(&vee));
            pt_v = pt_s;

            if nbits[i] {
                let pt_s = curve.add(&pt_v, pt_p);
                let ell = curve.line(&pt_v, pt_p, pt_q)?;
                let vee = curve.line(&pt_s, &curve.invert(&pt_s)?, pt_q)?;
                t = t.mul(&ell.div(&vee));
                pt_v = pt_s;
            }

            if i == 0 {
                break;
            }
            i -= 1;
        }
    }

    // Inversion for the Ate pairing
    if !sign {
        let vee = curve.line(&pt_v, &curve.invert(&pt_v)?, pt_q)?;
        t = t.mul(&vee).invert();
    }

    Ok(t)
}

/// Weil pairing
// /!\ I'm not checking that P, Q are on the same curve, I'm not checking that they are of the given order
// If you input incorrect data you get incorrect results
pub fn weil_pairing<F: Field + Clone + PartialEq>(
    curve: &EllipticCurve<F>,
    pt_p: ECPoint<F>,
    pt_q: ECPoint<F>,
    order: Integer,
) -> Result<F, ErrorKind> {
    let one = F::one();

    // P = Q, P = 0, or Q = 0
    if pt_p == pt_q || pt_p == ECPoint::PointAtInfinity || pt_q == ECPoint::PointAtInfinity {
        return Ok(one);
    }

    // Weil pairing
    let f_pq = miller(curve, &pt_p, &pt_q, &order)?;
    let f_qp = miller(curve, &pt_q, &pt_p, &order)?;
    let ratio = f_pq.div(&f_qp);

    // Sign correction if needed
    if order.is_odd() {
        Ok(ratio.neg())
    } else {
        Ok(ratio)
    }
}

/// Reduced Tate pairing
// /!\ I'm not checking that P is of the given order, or that the embedding degree is correct
// or that P and Q are on the same curve!
//
// Returns f_{n,P}(Q)^e  where div(f_{n,P}) = n(P) - n(O) and
// e = (q^k - 1)/n with q = base field size, n = order, and k = embedding degree
pub fn tate_pairing<F: Field + Clone + PartialEq>(
    curve: &EllipticCurve<F>,
    pt_p: &ECPoint<F>,
    pt_q: &ECPoint<F>,
    order: &Integer,
    embedding_degree: &Integer,
) -> Result<F, &'static str> {
    let q = F::base_order();

    // Check whether we need to move poles
    if let Ok(res) = miller(curve, pt_p, pt_q, order) {
        // We don't
        let one: Integer = 1.into();
        let e = q.large_pow(embedding_degree).sub(one).div(order);
        Ok(res.pow(&e))
    } else {
        // We do

        let pt_r = curve.clone().random_point();
        let f_qr = tate_pairing(
            curve,
            pt_p,
            &curve.add(pt_q, &pt_r),
            order,
            embedding_degree,
        )?;
        let f_r = tate_pairing(curve, pt_p, &pt_r, order, embedding_degree)?;

        Ok(f_qr.div(&f_r))
    }
}

/// N-th modified Ate pairing
// order is order of P and Q
// trace is the trace of the Frob over the base field
// trace_m_1 = trace - 1
// /!\ No checks
// P and Q on same curve
// P is on E/Fq
// Q is on E/Fq^k where k = embedding degree
// P in ker(Frob - 1)
// Q in ker(Frob - q)
pub fn ate_pairing<F: Field + Clone + PartialEq>(
    curve: &EllipticCurve<F>,
    pt_p: &ECPoint<F>,
    pt_q: &ECPoint<F>,
    order: &Integer,
    embedding_degree: &Integer,
    trace_m_1: &Integer,
) -> Result<F, ErrorKind> {
    let q = F::base_order();
    let res = miller(curve, pt_q, pt_p, trace_m_1)?;
    let one: Integer = 1.into();
    let e = q.large_pow(embedding_degree).sub(one).div(order);
    Ok(res.pow(&e))
}
