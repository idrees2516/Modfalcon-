use nalgebra::{DMatrix, DVector};
use rand_chacha::ChaCha20Rng;
use rand::SeedableRng;
use super::constants::*;
use super::ring::SecureRingElement;
use super::ntru::SecureNTRUBasis;
use super::sampler::SecureGaussianSampler;
use zeroize::ZeroizeOnDrop;

#[derive(ZeroizeOnDrop)]
pub struct SecureTrapdoorSampler {
    basis: SecureNTRUBasis,
    gs_basis: DMatrix<f64>,
    mu: DMatrix<f64>,
    sigma: f64,
    rng: ChaCha20Rng,
    sampler: SecureGaussianSampler,
}

impl SecureTrapdoorSampler {
    pub fn new(basis: SecureNTRUBasis) -> Self {
        let (gs_basis, mu) = Self::compute_gram_schmidt(&basis);
        
        SecureTrapdoorSampler {
            basis,
            gs_basis,
            mu,
            sigma: SIGMA,
            rng: ChaCha20Rng::from_entropy(),
            sampler: SecureGaussianSampler::new(SIGMA),
        }
    }

    fn compute_gram_schmidt(basis: &SecureNTRUBasis) -> (DMatrix<f64>, DMatrix<f64>) {
        let basis_matrix = basis.to_matrix();
        let (q, r) = basis_matrix.qr().unpack();
        let mu = &q.transpose() * &basis_matrix;
        (q, mu)
    }

    pub fn sample(&mut self, target: &SecureRingElement) -> Option<SecureRingElement> {
        let t = self.basis.canonical_embedding(target);
        let z = self.gaussian_lattice_sample(&t)?;
        
        let coeffs = z.iter()
            .map(|&x| x.round() as i64)
            .collect();
        
        Some(SecureRingElement::new(coeffs))
    }

    fn gaussian_lattice_sample(&mut self, target: &DVector<f64>) -> Option<DVector<f64>> {
        let n = self.mu.ncols();
        let mut result = DVector::zeros(n);
        
        for i in (0..n).rev() {
            let c = self.compute_center(target, &result, i);
            let sigma_i = self.compute_sigma(i);
            result[i] = self.sampler.sample_z() as f64 * sigma_i + c;
        }
        
        Some(result)
    }

    fn compute_center(&self, target: &DVector<f64>, partial: &DVector<f64>, i: usize) -> f64 {
        let mut center = target[i];
        for j in (i + 1)..self.mu.ncols() {
            center -= self.mu[(i, j)] * partial[j];
        }
        center /= self.mu[(i, i)];
        center
    }

    fn compute_sigma(&self, i: usize) -> f64 {
        let s = self.sigma;
        let r = (0..i).fold(0.0, |acc, j| {
            acc + (self.mu[(j, i)] / self.mu[(j, j)]).powi(2)
        });
        s * (1.0 - r).sqrt()
    }
}