use super::constants::*;

pub fn ntt_forward(coeffs: &mut [i64]) {
    let n = coeffs.len();
    let mut m = 1;
    while m < n {
        let h = n / (2 * m);
        for i in 0..m {
            let w = NTT_ROOTS[m + i];
            for j in 0..h {
                let k = 2 * j * m + i;
                let u = coeffs[k];
                let v = (coeffs[k + m] * w) % Q;
                coeffs[k] = (u + v) % Q;
                coeffs[k + m] = (u - v + Q) % Q;
            }
        }
        m *= 2;
    }
}

pub fn ntt_inverse(coeffs: &mut [i64]) {
    let n = coeffs.len();
    let mut m = n / 2;
    while m >= 1 {
        let h = n / (2 * m);
        for i in 0..m {
            let w = NTT_ROOTS_INV[m + i];
            for j in 0..h {
                let k = 2 * j * m + i;
                let u = coeffs[k];
                let v = coeffs[k + m];
                coeffs[k] = (u + v) % Q;
                coeffs[k + m] = ((u - v) * w) % Q;
            }
        }
        m /= 2;
    }

    let n_inv = mod_inverse(n as i64, Q);
    for coeff in coeffs.iter_mut() {
        *coeff = (*coeff * n_inv) % Q;
    }
}

fn mod_inverse(a: i64, m: i64) -> i64 {
    let mut t = 0i64;
    let mut newt = 1i64;
    let mut r = m;
    let mut newr = a;

    while newr != 0 {
        let quotient = r / newr;
        (t, newt) = (newt, t - quotient * newt);
        (r, newr) = (newr, r - quotient * newr);
    }

    if r > 1 {
        panic!("a is not invertible");
    }
    if t < 0 {
        t += m;
    }
    t
}