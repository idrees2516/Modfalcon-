use rand::RngCore;
use rand_chacha::ChaCha20Rng;
use rand::SeedableRng;
use super::constants::*;
use super::ring::SecureRingElement;
use zeroize::ZeroizeOnDrop;

#[derive(ZeroizeOnDrop)]
pub struct SecureGaussianSampler {
    sigma: f64,
    precision_bits: usize,
    cdt: Vec<u64>,
    rng: ChaCha20Rng,
}

impl SecureGaussianSampler {
    pub fn new(sigma: f64) -> Self {
        let precision_bits = 128;
        let cdt = Self::precompute_cdt(sigma, precision_bits);
        
        SecureGaussianSampler {
            sigma,
            precision_bits,
            cdt,
            rng: ChaCha20Rng::from_entropy(),
        }
    }

    fn precompute_cdt(sigma: f64, precision_bits: usize) -> Vec<u64> {
        let tail = (12.0 * sigma).ceil() as i64;
        let mut cdt = Vec::with_capacity((2 * tail + 1) as usize);
        let scale = (1u128 << precision_bits) as f64;
        
        let mut sum = 0.0;
        for x in -tail..=tail {
            let prob = (-((x as f64).powi(2)) / (2.0 * sigma.powi(2))).exp();
            sum += prob;
            cdt.push((sum * scale) as u64);
        }
        
        cdt
    }

    pub fn sample_z(&mut self) -> i64 {
        let u = self.rng.next_u64();
        let mut result = 0i64;
        
        for (i, &threshold) in self.cdt.iter().enumerate() {
            let is_less = (u <= threshold) as i64;
            result = (result & !is_less) | (i as i64 & is_less);
        }
        
        result - (self.cdt.len() as i64 / 2)
    }

    pub fn sample_ring(&mut self) -> SecureRingElement {
        let coeffs: Vec<i64> = (0..D)
            .map(|_| self.sample_z())
            .collect();
        SecureRingElement::new(coeffs)
    }
}