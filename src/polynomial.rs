use super::constants::Q;

#[derive(Clone, Debug)]
pub struct Polynomial {
    pub coeffs: Vec<i128>,
}

impl Polynomial {
    pub fn from_coeffs(coeffs: &[i128]) -> Self {
        Polynomial {
            coeffs: coeffs.to_vec(),
        }
    }

    pub fn degree(&self) -> usize {
        self.coeffs.len().saturating_sub(1)
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs.iter().all(|&c| c == 0)
    }

    pub fn add(&self, other: &Self, q: i128) -> Self {
        let max_len = self.coeffs.len().max(other.coeffs.len());
        let mut result = vec![0i128; max_len];
        
        for i in 0..max_len {
            let a = if i < self.coeffs.len() { self.coeffs[i] } else { 0 };
            let b = if i < other.coeffs.len() { other.coeffs[i] } else { 0 };
            result[i] = (a + b).rem_euclid(q);
        }
        
        Polynomial::from_coeffs(&result)
    }

    pub fn sub(&self, other: &Self, q: i128) -> Self {
        let max_len = self.coeffs.len().max(other.coeffs.len());
        let mut result = vec![0i128; max_len];
        
        for i in 0..max_len {
            let a = if i < self.coeffs.len() { self.coeffs[i] } else { 0 };
            let b = if i < other.coeffs.len() { other.coeffs[i] } else { 0 };
            result[i] = (a - b).rem_euclid(q);
        }
        
        Polynomial::from_coeffs(&result)
    }

    pub fn mul(&self, other: &Self, q: i128) -> Self {
        let n = self.coeffs.len() + other.coeffs.len() - 1;
        let mut result = vec![0i128; n];
        
        for (i, &a) in self.coeffs.iter().enumerate() {
            for (j, &b) in other.coeffs.iter().enumerate() {
                result[i + j] = (result[i + j] + a * b).rem_euclid(q);
            }
        }
        
        Polynomial::from_coeffs(&result)
    }

    pub fn div_rem(&self, other: &Self, q: i128) -> Option<(Self, Self)> {
        if other.is_zero() {
            return None;
        }

        let mut quotient = vec![0i128; self.degree().saturating_sub(other.degree()) + 1];
        let mut remainder = self.coeffs.clone();
        
        while remainder.len() >= other.coeffs.len() {
            let lead_div = mod_inverse(other.coeffs[other.coeffs.len()-1], q)?;
            let factor = (remainder[remainder.len()-1] * lead_div).rem_euclid(q);
            let pos = remainder.len() - other.coeffs.len();
            quotient[pos] = factor;
            
            for (i, &b) in other.coeffs.iter().enumerate() {
                remainder[pos + i] = (remainder[pos + i] - factor * b).rem_euclid(q);
            }
            
            while !remainder.is_empty() && remainder[remainder.len()-1] == 0 {
                remainder.pop();
            }
        }
        
        Some((
            Polynomial::from_coeffs(&quotient),
            Polynomial::from_coeffs(&remainder)
        ))
    }

    pub fn mod_inv(&self, modulus: &Self, q: i128) -> Option<Self> {
        let (mut s_prev, mut s) = (Polynomial::from_coeffs(&[1]), Polynomial::from_coeffs(&[0]));
        let (mut t_prev, mut t) = (Polynomial::from_coeffs(&[0]), Polynomial::from_coeffs(&[1]));
        let (mut r_prev, mut r) = (self.clone(), modulus.clone());
        
        while !r.is_zero() {
            let (quotient, remainder) = r_prev.div_rem(&r, q)?;
            r_prev = r;
            r = remainder;
            
            let s_new = s_prev.sub(&quotient.mul(&s, q), q);
            s_prev = s;
            s = s_new;
            
            let t_new = t_prev.sub(&quotient.mul(&t, q), q);
            t_prev = t;
            t = t_new;
        }
        
        if r_prev.degree() == 0 {
            let scale = mod_inverse(r_prev.coeffs[0], q)?;
            Some(s_prev.mul(&Polynomial::from_coeffs(&[scale]), q))
        } else {
            None
        }
    }
}

fn mod_inverse(a: i128, m: i128) -> Option<i128> {
    let mut t = 0i128;
    let mut newt = 1i128;
    let mut r = m;
    let mut newr = a;

    while newr != 0 {
        let quotient = r / newr;
        (t, newt) = (newt, t - quotient * newt);
        (r, newr) = (newr, r - quotient * newr);
    }

    if r > 1 {
        return None;
    }
    Some(if t < 0 { t + m } else { t })
}