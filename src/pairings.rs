use crate::{bigint::BigInt, elliptic_curve::ECPoint, field::Field};

// Miller's algorithm
// Returns f_{n,P}(Q) where div(f_{n,P}) = n(P) - ([n]P) - (n-1)(0)
pub fn miller<F: Field + Clone + Eq>(
    pt_p: &ECPoint<F>,
    pt_q: &ECPoint<F>,
    n: &BigInt,
) -> Result<F, &'static str> {
    // Basic checks
    if pt_p.is_zero() {
        return Err("P must not be zero");
    }
    if pt_q.is_zero() {
        return Err("Q must not be zero");
    }
    if n.is_zero() {
        return Ok(F::one());
    }

    // Negative values of n are allowed, in which case
    // Q is evaluated instead at (v_{[n]P} f_{n,P)})^(-1)
    let sign = n.is_positive();
    let n = n.abs();
    let nbits = n.to_bits();

    let one = F::one();

    let mut t = one;
    let mut i: usize = nbits.len() - 1; // Will not underflow because n != 0
    let mut pt_v = pt_p.clone();

    if i != 0 {
        i -= 1;

        // Miller loop
        loop {
            let pt_s = pt_v.double();
            let ell = pt_v.line(&pt_v, &pt_q);
            let vee = pt_s.line(&pt_s.invert(), &pt_q);
            t = t.square().mul(&ell.div(&vee));
            pt_v = pt_s;

            if nbits[i] {
                let pt_s = pt_v.add(&pt_p);
                let ell = pt_v.line(&pt_p, &pt_q);
                let vee = pt_s.line(&pt_s.invert(), &pt_q);
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
        let vee = pt_v.line(&pt_v.invert(), &pt_q);
        t = t.mul(&vee).invert();
    }

    Ok(t)
}

// Weil pairing
// /!\ I'm not checking that P, Q are on the same curve, I'm not checking that they are of the given order
// If you input incorrect data you get incorrect results
pub fn weil_pairing<F: Field + Clone + Eq>(
    pt_p: ECPoint<F>,
    pt_q: ECPoint<F>,
    order: BigInt,
) -> Result<F, &'static str> {
    let one = F::one();

    // P = Q, P = 0, or Q = 0
    if pt_p.is_equal(&pt_q) || pt_p.is_zero() || pt_q.is_zero() {
        return Ok(one);
    }

    // Weil pairing
    let f_pq = miller(&pt_p, &pt_q, &order)?;
    let f_qp = miller(&pt_q, &pt_p, &order)?;
    let ratio = f_pq.div(&f_qp);

    // Sign correction if needed
    if order.is_odd() {
        return Ok(ratio.neg());
    } else {
        return Ok(ratio);
    }
}

// Reduced Tate pairing
// /!\ I'm not checking that P is of the given order, or that the embedding degree is correct
// or that P and Q are on the same curve!
//
// Returns f_{n,P}(Q)^e  where div(f_{n,P}) = n(P) - n(O) and
// e = (q^k - 1)/n with q = base field size, n = order, and k = embedding degree
pub fn tate_pairing<F: Field + Clone + Eq>(
    pt_p: &ECPoint<F>,
    pt_q: &ECPoint<F>,
    order: &BigInt,
    embedding_degree: &BigInt,
) -> Result<F, &'static str> {
    let q = F::base_order();

    // Check whether we need to move poles
    if let Ok(res) = miller(&pt_p, &pt_q, &order) {
        // We don't
        let e = q.pow(embedding_degree).sub(&BigInt::one()).div(order);
        Ok(res.pow(&e))
    } else {
        // We do
        let pt_r = pt_p.clone().curve.random_point();
        let f_qr = tate_pairing(&pt_p, &pt_q.add(&pt_r), order, embedding_degree)?;
        let f_r = tate_pairing(&pt_p, &pt_r, order, embedding_degree)?;

        Ok(f_qr.div(&f_r))
    }
}

// N-th modified Ate pairing
// order is order of P and Q
// trace is the trace of the Frob over the base field
// trace_m_1 = trace - 1
// /!\ No checks
// P and Q on same curve
// P is on E/Fq
// Q is on E/Fq^k where k = embedding degree
// P in ker(Frob - 1)
// Q in ker(Frob - q)
pub fn ate_pairing<F: Field + Clone + Eq>(
    pt_p: &ECPoint<F>,
    pt_q: &ECPoint<F>,
    order: &BigInt,
    embedding_degree: &BigInt,
    trace_m_1: &BigInt,
) -> Result<F, &'static str> {
    let q = F::base_order();
    let res = miller(&pt_q, &pt_p, &trace_m_1)?;
    let e = q.pow(&embedding_degree).sub(&BigInt::one()).div(&order);
    Ok(res.pow(&e))
}
