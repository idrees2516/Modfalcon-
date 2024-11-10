use super::constants::*;
use super::ring::SecureRingElement;
use super::sampler::SecureGaussianSampler;
use nalgebra::{DMatrix, DVector};
use num_complex::Complex;
use rayon::prelude::*;
use std::f64::consts::PI;
use zeroize::ZeroizeOnDrop;

#[derive(Clone, ZeroizeOnDrop)]
pub struct SecureNTRUBasis {
    pub f: SecureRingElement,
    pub g: SecureRingElement,
    pub F: SecureRingElement,
    pub G: SecureRingElement,
    pub h: SecureRingElement,
    gram_schmidt_cache: Option<(DMatrix<f64>, Vec<f64>)>,
}

impl SecureNTRUBasis {
    pub fn generate() -> Option<Self> {
        let mut sampler = SecureGaussianSampler::new(SIGMA);
        let mut best_basis = None;
        let mut min_gs_norm = f64::INFINITY;

        for _ in 0..MAX_ITER {
            let f = sampler.sample_ring();
            let g = sampler.sample_ring();

            if let Some(f_inv_q) = f.invert() {
                let h = f_inv_q.multiply(&g);
                let (F, G) = Self::compute_short_vectors(&f, &g)?;
                
                let basis = SecureNTRUBasis {
                    f: f.clone(),
                    g: g.clone(),
                    F,
                    G,
                    h,
                    gram_schmidt_cache: None,
                };

                if basis.verify_basis() {
                    let gs_norm = basis.compute_max_gs_norm();
                    if gs_norm < min_gs_norm {
                        min_gs_norm = gs_norm;
                        best_basis = Some(basis);
                    }
                }
            }
        }
        best_basis
    }

    fn compute_short_vectors(f: &SecureRingElement, g: &SecureRingElement) -> Option<(SecureRingElement, SecureRingElement)> {
        let mut F_coeffs = vec![0i64; D];
        let mut G_coeffs = vec![0i64; D];
        
        // Extended GCD over polynomial ring
        let (gcd, s, t) = Self::xgcd_poly(f, g)?;
        
        // Scale to achieve fG - gF = q
        let scale = Q / gcd.coeffs[0];
        F_coeffs.copy_from_slice(&s.coeffs);
        G_coeffs.copy_from_slice(&t.coeffs);
        
        for coeff in F_coeffs.iter_mut() {
            *coeff = (*coeff * scale).rem_euclid(Q);
        }
        for coeff in G_coeffs.iter_mut() {
            *coeff = (*coeff * scale).rem_euclid(Q);
        }
        
        Some((
            SecureRingElement::new(F_coeffs),
            SecureRingElement::new(G_coeffs)
        ))
    }

    fn xgcd_poly(f: &SecureRingElement, g: &SecureRingElement) -> Option<(SecureRingElement, SecureRingElement, SecureRingElement)> {
        let mut r0 = f.clone();
        let mut r1 = g.clone();
        let mut s0 = SecureRingElement::one();
        let mut s1 = SecureRingElement::zero();
        let mut t0 = SecureRingElement::zero();
        let mut t1 = SecureRingElement::one();

        while !r1.coeffs.iter().all(|&x| x == 0) {
            let (q, r) = Self::poly_divmod(&r0, &r1)?;
            r0 = r1;
            r1 = r;

            let new_s = s0.sub(&q.multiply(&s1));
            s0 = s1;
            s1 = new_s;

            let new_t = t0.sub(&q.multiply(&t1));
            t0 = t1;
            t1 = new_t;
        }

        Some((r0, s0, t0))
    }

    fn poly_divmod(a: &SecureRingElement, b: &SecureRingElement) -> Option<(SecureRingElement, SecureRingElement)> {
        if b.coeffs.iter().all(|&x| x == 0) {
            return None;
        }

        let mut q_coeffs = vec![0i64; D];
        let mut r = a.clone();
        let b_lead = b.coeffs.last().copied()?;
        let b_deg = b.coeffs.len() - 1;

        while !r.coeffs.is_empty() && r.coeffs.len() > b_deg {
            let r_lead = *r.coeffs.last()?;
            let q_coeff = (r_lead * Self::mod_inverse(b_lead, Q)?).rem_euclid(Q);
            let pos = r.coeffs.len() - b_deg - 1;
            q_coeffs[pos] = q_coeff;

            for (i, &b_coeff) in b.coeffs.iter().enumerate() {
                r.coeffs[pos + i] = (r.coeffs[pos + i] - q_coeff * b_coeff).rem_euclid(Q);
            }

            while !r.coeffs.is_empty() && r.coeffs[r.coeffs.len() - 1] == 0 {
                r.coeffs.pop();
            }
        }

        Some((SecureRingElement::new(q_coeffs), r))
    }

    fn mod_inverse(a: i64, m: i64) -> Option<i64> {
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
            return None;
        }
        Some(if t < 0 { t + m } else { t })
    }

    pub fn verify_basis(&self) -> bool {
        // Verify fG - gF = q
        let fG = self.f.multiply(&self.G);
        let gF = self.g.multiply(&self.F);
        let det = fG.sub(&gF);
        
        let mut q_poly = SecureRingElement::zero();
        q_poly.coeffs[0] = Q;
        
        if det != q_poly {
            return false;
        }

        // Verify Gram-Schmidt norms
        let gs_norms = self.compute_gram_schmidt_norms();
        let bound = RHO * SIGMA * (2.0 * D as f64).sqrt();

        gs_norms.iter().all(|&norm| norm <= bound)
    }

    pub fn compute_gram_schmidt_norms(&self) -> Vec<f64> {
        if let Some((_, norms)) = &self.gram_schmidt_cache {
            return norms.clone();
        }

        let basis_matrix = self.to_matrix();
        let (q, _) = basis_matrix.qr().unpack();
        let norms: Vec<f64> = q.column_iter()
            .map(|col| col.norm())
            .collect();

        norms
    }

    pub fn compute_max_gs_norm(&self) -> f64 {
        self.compute_gram_schmidt_norms()
            .into_iter()
            .fold(0.0, f64::max)
    }

    pub fn to_matrix(&self) -> DMatrix<f64> {
        let n = 2 * D;
        let mut matrix = DMatrix::zeros(n, n);

        // Embed f and g
        for i in 0..D {
            let angle = 2.0 * PI * i as f64 / D as f64;
            let root = Complex::new(angle.cos(), angle.sin());
            
            let f_val = self.evaluate_at(&self.f, root);
            let g_val = self.evaluate_at(&self.g, root);
            
            matrix[(2*i, 0)] = f_val.re;
            matrix[(2*i+1, 0)] = f_val.im;
            matrix[(2*i, 1)] = g_val.re;
            matrix[(2*i+1, 1)] = g_val.im;
        }

        // Embed F and G
        for i in 0..D {
            let angle = 2.0 * PI * i as f64 / D as f64;
            let root = Complex::new(angle.cos(), angle.sin());
            
            let F_val = self.evaluate_at(&self.F, root);
            let G_val = self.evaluate_at(&self.G, root);
            
            matrix[(2*i, 2)] = F_val.re;
            matrix[(2*i+1, 2)] = F_val.im;
            matrix[(2*i, 3)] = G_val.re;
            matrix[(2*i+1, 3)] = G_val.im;
        }

        matrix
    }

    fn evaluate_at(&self, poly: &SecureRingElement, z: Complex<f64>) -> Complex<f64> {
        let mut result = Complex::new(0.0, 0.0);
        let mut power = Complex::new(1.0, 0.0);
        
        for &coeff in &poly.coeffs {
            result += power * (coeff as f64);
            power *= z;
        }
        
        result
    }

    pub fn canonical_embedding(&self, element: &SecureRingElement) -> DVector<f64> {
        let mut embedding = DVector::zeros(2 * D);
        
        for i in 0..D {
            let angle = 2.0 * PI * i as f64 / D as f64;
            let root = Complex::new(angle.cos(), angle.sin());
            let value = self.evaluate_at(element, root);
            
            embedding[2*i] = value.re;
            embedding[2*i+1] = value.im;
        }
        
        embedding
    }
}