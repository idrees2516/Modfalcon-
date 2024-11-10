use super::constants::*;
use super::ring::SecureRingElement;
use super::ntru::SecureNTRUBasis;
use super::sampler::SecureGaussianSampler;
use super::trapdoor::SecureTrapdoorSampler;
use sha3::{Digest, Sha3_512};
use subtle::{Choice, ConstantTimeEq};
use zeroize::{Zeroize, ZeroizeOnDrop};
use rand_chacha::ChaCha20Rng;
use rand::SeedableRng;
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(ZeroizeOnDrop)]
pub struct SecureModFalcon {
    basis: SecureNTRUBasis,
    sampler: SecureTrapdoorSampler,
    salt_counter: AtomicU64,
    rng: ChaCha20Rng,
}

#[derive(Clone, Debug, ZeroizeOnDrop)]
pub struct Signature {
    s1: SecureRingElement,
    s2: SecureRingElement,
    salt: [u8; 32],
    nonce: [u8; 32],
}

impl SecureModFalcon {
    pub fn new() -> Option<Self> {
        let basis = SecureNTRUBasis::generate()?;
        let sampler = SecureTrapdoorSampler::new(basis.clone());
        
        Some(SecureModFalcon {
            basis,
            sampler,
            salt_counter: AtomicU64::new(0),
            rng: ChaCha20Rng::from_entropy(),
        })
    }

    pub fn sign(&mut self, message: &[u8]) -> Option<Signature> {
        let mut salt = [0u8; 32];
        let mut nonce = [0u8; 32];
        self.rng.fill_bytes(&mut salt);
        self.rng.fill_bytes(&mut nonce);

        let counter = self.salt_counter.fetch_add(1, Ordering::SeqCst);
        let hashed_message = self.hash_to_ring(message, &salt, &nonce, counter);
        
        let mut attempts = 0;
        while attempts < MAX_ITER {
            if let Some(s1) = self.sampler.sample(&hashed_message) {
                let s2 = self.compute_s2(&s1, &hashed_message)?;
                if self.verify_norm(&s1, &s2) {
                    return Some(Signature { s1, s2, salt, nonce });
                }
            }
            attempts += 1;
        }
        None
    }

    pub fn verify(&self, message: &[u8], signature: &Signature) -> bool {
        let counter = self.salt_counter.load(Ordering::SeqCst);
        let hashed_message = self.hash_to_ring(message, &signature.salt, &signature.nonce, counter);
        
        // Verify norm bounds
        if !self.verify_norm(&signature.s1, &signature.s2) {
            return false;
        }

        // Verify s1 + s2*h = c mod q
        let s2h = signature.s2.multiply(&self.basis.h);
        let lhs = signature.s1.add(&s2h);
        lhs.ct_eq(&hashed_message).unwrap_u8() == 1
    }

    fn compute_s2(&self, s1: &SecureRingElement, c: &SecureRingElement) -> Option<SecureRingElement> {
        let diff = c.sub(s1);
        let h_inv = self.basis.h.invert()?;
        Some(diff.multiply(&h_inv))
    }

    fn verify_norm(&self, s1: &SecureRingElement, s2: &SecureRingElement) -> bool {
        let norm_bound = BETA as f64 * SIGMA;
        s1.norm() <= norm_bound && s2.norm() <= norm_bound
    }

    fn hash_to_ring(&self, message: &[u8], salt: &[u8; 32], nonce: &[u8; 32], counter: u64) -> SecureRingElement {
        let mut hasher = Sha3_512::new();
        hasher.update(salt);
        hasher.update(nonce);
        hasher.update(counter.to_be_bytes());
        hasher.update(message);
        let hash = hasher.finalize();

        let mut coeffs = Vec::with_capacity(D);
        let mut block_counter = 0u64;
        
        while coeffs.len() < D {
            let mut block_hasher = Sha3_512::new();
            block_hasher.update(&hash);
            block_hasher.update(block_counter.to_be_bytes());
            let block = block_hasher.finalize();
            
            for chunk in block.chunks(8) {
                if coeffs.len() >= D {
                    break;
                }
                let mut value = 0i64;
                for &byte in chunk {
                    value = (value << 8) | (byte as i64);
                }
                coeffs.push(value.rem_euclid(Q));
            }
            block_counter += 1;
        }

        SecureRingElement::new(coeffs)
    }

    pub fn public_key(&self) -> SecureRingElement {
        self.basis.h.clone()
    }

    pub fn verify_key_pair(&self) -> bool {
        self.basis.verify_basis()
    }
}

impl Zeroize for Signature {
    fn zeroize(&mut self) {
        self.s1.zeroize();
        self.s2.zeroize();
        self.salt.zeroize();
        self.nonce.zeroize();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sign_verify() {
        let mut falcon = SecureModFalcon::new().unwrap();
        let message = b"test message";
        let signature = falcon.sign(message).unwrap();
        assert!(falcon.verify(message, &signature));
    }

    #[test]
    fn test_signature_rejection() {
        let mut falcon = SecureModFalcon::new().unwrap();
        let message = b"test message";
        let signature = falcon.sign(message).unwrap();
        let wrong_message = b"wrong message";
        assert!(!falcon.verify(wrong_message, &signature));
    }

    #[test]
    fn test_key_pair_verification() {
        let falcon = SecureModFalcon::new().unwrap();
        assert!(falcon.verify_key_pair());
    }

    #[test]
    fn test_norm_bounds() {
        let mut falcon = SecureModFalcon::new().unwrap();
        let message = b"test message";
        let signature = falcon.sign(message).unwrap();
        assert!(falcon.verify_norm(&signature.s1, &signature.s2));
    }
}