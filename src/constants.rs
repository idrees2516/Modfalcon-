use std::f64::consts::PI;

pub const SECURITY_LEVEL: usize = 128;
pub const D: usize = 512;
pub const Q: i64 = 12289;
pub const RHO: f64 = 1.05;
pub const MAX_ITER: usize = 10000;
pub const SIGMA: f64 = 1.17 * (2.0 * PI).sqrt();
pub const BETA: usize = 380;
pub const REJECTION_BOUND: f64 = 1.47;

pub const NTT_ROOTS: [i64; D] = compute_ntt_roots();
pub const NTT_ROOTS_INV: [i64; D] = compute_ntt_roots_inv();

const fn compute_ntt_roots() -> [i64; D] {
    let mut roots = [0i64; D];
    let primitive_root = find_primitive_root(Q);
    let factor = mod_pow(primitive_root, (Q - 1) / (2 * D) as i64, Q);
    
    let mut current = 1;
    let mut i = 0;
    while i < D {
        roots[i] = current;
        current = (current * factor) % Q;
        i += 1;
    }
    roots
}

const fn compute_ntt_roots_inv() -> [i64; D] {
    let mut roots = [0i64; D];
    let primitive_root = find_primitive_root(Q);
    let factor = mod_pow(primitive_root, (Q - 1) - (Q - 1) / (2 * D) as i64, Q);
    
    let mut current = 1;
    let mut i = 0;
    while i < D {
        roots[i] = current;
        current = (current * factor) % Q;
        i += 1;
    }
    roots
}

const fn find_primitive_root(q: i64) -> i64 {
    let mut g = 2;
    while !is_primitive_root(g, q) {
        g += 1;
    }
    g
}

const fn is_primitive_root(g: i64, q: i64) -> bool {
    let factors = [2i64, (q-1)/2];
    let mut i = 0;
    while i < factors.len() {
        if mod_pow(g, (q-1)/factors[i], q) == 1 {
            return false;
        }
        i += 1;
    }
    true
}

const fn mod_pow(mut base: i64, mut exp: i64, modulus: i64) -> i64 {
    let mut result = 1;
    base = base % modulus;
    while exp > 0 {
        if exp & 1 == 1 {
            result = (result * base) % modulus;
        }
        base = (base * base) % modulus;
        exp >>= 1;
    }
    result
}