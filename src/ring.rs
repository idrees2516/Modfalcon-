use super::constants::*;
use super::polynomial::Polynomial;
use rayon::prelude::*;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};
use zeroize::ZeroizeOnDrop;

#[derive(Clone, Debug, PartialEq, ZeroizeOnDrop)]
pub struct SecureRingElement {
    coeffs: Vec<i64>,
}

impl SecureRingElement {
    pub fn new(coeffs: Vec<i64>) -> Self {
        assert_eq!(coeffs.len(), D);
        let reduced_coeffs = coeffs.into_iter()
            .map(|x| x.rem_euclid(Q))
            .collect();
        SecureRingElement { coeffs: reduced_coeffs }
    }

    pub fn zero() -> Self {
        SecureRingElement { coeffs: vec![0; D] }
    }

    pub fn one() -> Self {
        let mut coeffs = vec![0; D];
        coeffs[0] = 1;
        SecureRingElement::new(coeffs)
    }

    pub fn random(rng: &mut impl rand::RngCore) -> Self {
        let coeffs: Vec<i64> = (0..D)
            .map(|_| rng.next_u64() as i64 % Q)
            .collect();
        SecureRingElement::new(coeffs)
    }

    pub fn add(&self, other: &Self) -> Self {
        let coeffs = self.coeffs.par_iter()
            .zip(&other.coeffs)
            .map(|(&a, &b)| (a + b).rem_euclid(Q))
            .collect();
        SecureRingElement::new(coeffs)
    }

    pub fn sub(&self, other: &Self) -> Self {
        let coeffs = self.coeffs.par_iter()
            .zip(&other.coeffs)
            .map(|(&a, &b)| (a - b).rem_euclid(Q))
            .collect();
        SecureRingElement::new(coeffs)
    }

    pub fn multiply(&self, other: &Self) -> Self {
        let mut result = vec![0i64; 2 * D];
        
        self.coeffs.par_iter().enumerate().for_each(|(i, &a)| {
            other.coeffs.iter().enumerate().for_each(|(j, &b)| {
                let idx = (i + j) % D;
                let mut temp = result[idx];
                temp = (temp + a * b).rem_euclid(Q);
                result[idx] = temp;
            });
        });

        for i in D..(2 * D) {
            result[i - D] = (result[i - D] - result[i]).rem_euclid(Q);
        }
        result.truncate(D);
        
        SecureRingElement::new(result)
    }

    pub fn scalar_mult(&self, scalar: i64) -> Self {
        let coeffs = self.coeffs.par_iter()
            .map(|&x| (x * scalar).rem_euclid(Q))
            .collect();
        SecureRingElement::new(coeffs)
    }

    pub fn invert(&self) -> Option<Self> {
        let modulus = Polynomial::from_coeffs(&[1] + &vec![0; D-1] + &[1]);
        let poly = Polynomial::from_coeffs(&self.coeffs.iter().map(|&x| x as i128).collect::<Vec<_>>());
        
        let inv_poly = poly.mod_inv(&modulus, Q as i128)?;
        let coeffs = inv_poly.coeffs.iter()
            .map(|&x| x.rem_euclid(Q as i128) as i64)
            .collect();
        
        Some(SecureRingElement::new(coeffs))
    }

    pub fn norm_squared(&self) -> f64 {
        self.coeffs.par_iter()
            .map(|&x| {
                let reduced = if x > Q/2 { x - Q } else { x };
                (reduced as f64).powi(2)
            })
            .sum()
    }

    pub fn norm(&self) -> f64 {
        self.norm_squared().sqrt()
    }

    pub fn to_ntt(&self) -> Self {
        let mut coeffs = self.coeffs.clone();
        ntt_forward(&mut coeffs);
        SecureRingElement::new(coeffs)
    }

    pub fn from_ntt(&self) -> Self {
        let mut coeffs = self.coeffs.clone();
        ntt_inverse(&mut coeffs);
        SecureRingElement::new(coeffs)
    }
}

impl ConditionallySelectable for SecureRingElement {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        let coeffs = a.coeffs.iter().zip(&b.coeffs)
            .map(|(&x, &y)| {
                let x_int = x as i128;
                let y_int = y as i128;
                let selected = i128::conditional_select(&x_int, &y_int, choice);
                selected as i64
            })
            .collect();
        SecureRingElement::new(coeffs)
    }
}

impl ConstantTimeEq for SecureRingElement {
    fn ct_eq(&self, other: &Self) -> Choice {
        let mut result = Choice::from(1u8);
        for (a, b) in self.coeffs.iter().zip(&other.coeffs) {
            result &= Choice::from((a == b) as u8);
        }
        result
    }
}