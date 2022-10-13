use std::ops::MulAssign;
use ark_ec::ProjectiveCurve;
use ark_ff::{One, PrimeField, UniformRand};
use ark_std::rand::Rng;


use crate::{g1_base, g2_base, Scalar, G1, G2};

/// Structured Reference String, defined over BLS12-381 curve.
#[derive(Debug, Clone)]
pub struct SRS {
    pub g1: Vec<G1>,
    pub g2: Vec<G2>,
}

impl SRS {
    pub fn new<R: Rng>(max_degree: usize, rng: &mut R) -> SRS {
        let r = Scalar::rand(rng);

        let mut g1 = Vec::new();
        g1.push(g1_base());
        for i in 0..max_degree {
            let ele = g1[i].mul(r.into_repr());
            g1.push(ele)
        }

        let g2 = vec![g2_base(), g2_base().mul(r.into_repr())];

        SRS { g1, g2 }
    }

    /// Update SRS.
    pub fn update<R: Rng>(&mut self, rng: &mut R) {
        assert!(self.g2.len() == 2);
        let r = Scalar::rand(rng);

        let mut r_pow = Scalar::one();
        for x in self.g1.iter_mut().skip(1) {
            r_pow.mul_assign(&r);
            x.mul_assign(r_pow);
        }

        self.g2[1].mul_assign(r);
    }

    /// The public setup parameters come from https://github.com/findora-crypto/export-setup-parameters.
    pub fn load_from_public_setup_parameters(max_degree: usize) -> SRS {
        let g1 = pub_srs::export_g1_from_public_setup_parameters(max_degree);
        let g2 = pub_srs::export_g2_from_public_setup_parameters();

        SRS { g1, g2 }
    }
}
