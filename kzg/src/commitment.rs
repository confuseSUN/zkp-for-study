use ark_bls12_381::Bls12_381;
use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_ff::{One, PrimeField};
use ark_poly::{univariate::DenseOrSparsePolynomial, Polynomial, UVPolynomial};
use std::ops::{Neg, Sub};

use crate::{srs::SRS, Scalar, G1};

type Poly = ark_poly::polynomial::univariate::DensePolynomial<Scalar>;

pub struct KZGCommitmentScheme<'a>(&'a SRS);

pub struct KZGCommitmentProof {
    opening_value: Scalar,
    comm_h: G1,
}

impl KZGCommitmentScheme<'_> {
    fn max_degree(&self) -> usize {
        self.0.g1.len() - 1
    }

    pub fn commit(&self, poly: &Poly) -> G1 {
        assert!(self.max_degree() >= poly.degree());

        let g1 = &self.0.g1;
        let coefs = &poly.coeffs;

        // The msm(Multi-scalar multiplication) fast algorithm should be used.
        let commitment = coefs
            .iter()
            .zip(g1.iter())
            .map(|(coef, g)| g.mul(coef.into_repr()))
            .sum::<G1>();
        commitment
    }

    pub fn prove(&self, poly: &Poly, z: Scalar) -> KZGCommitmentProof {
        let opening_value = poly.evaluate(&z);

        let dividend = poly.sub(&Poly::from_coefficients_slice(&[opening_value]));
        let divisor = Poly::from_coefficients_slice(&[z.neg(), Scalar::one()]);

        let (quotient, remainder) =
            DenseOrSparsePolynomial::divide_with_q_and_r(&(&dividend).into(), &(&divisor).into())
                .unwrap();

        assert!(remainder.is_empty());

        let comm_h = self.commit(&quotient);

        KZGCommitmentProof {
            opening_value,
            comm_h,
        }
    }

    pub fn verify(&self, poly_comm: &G1, proof: &KZGCommitmentProof, z: &Scalar) -> bool {
        let g1_base = self.0.g1[0];
        let g2_base = self.0.g2[0];
        let g2_r = self.0.g2[1];

        let left = Bls12_381::pairing(
            poly_comm.sub(g1_base.mul(proof.opening_value.into_repr())),
            g2_base,
        );
        let right = Bls12_381::pairing(proof.comm_h, g2_r.sub(g2_base.mul(z.into_repr())));

        left == right
    }
}

#[cfg(test)]
mod test_kzg {
    use super::{KZGCommitmentScheme, Poly, SRS};
    use crate::Scalar;
    use ark_ff::UniformRand;
    use ark_poly::UVPolynomial;
    use ark_std::{rand::Rng, test_rng};

    #[test]
    fn test_kzg_comm() {
        let max_degree = 20;
        let mut rng = test_rng();
        let mut srs = SRS::new(max_degree, &mut rng);
        kzg_comm(srs.clone(), max_degree, &mut rng);

        srs.update(&mut rng);
        kzg_comm(srs, max_degree, &mut rng);
    }

    #[ignore = "git lfs disable"]
    #[test]
    fn test_kzg_comm_with_public_srs() {
        let max_degree = 20;
        let mut rng = test_rng();
        let mut pub_srs = SRS::load_from_public_setup_parameters(max_degree);
        kzg_comm(pub_srs.clone(), max_degree, &mut rng);

        pub_srs.update(&mut rng);
        kzg_comm(pub_srs, max_degree, &mut rng);
    }

    fn kzg_comm<R: Rng>(srs: SRS, max_degree: usize, rng: &mut R) {
        let kzg_comm_scheme = KZGCommitmentScheme(&srs);

        let mut coefs = Vec::new();
        for _ in 0..max_degree + 1 {
            let coef = Scalar::rand(rng);
            coefs.push(coef);
        }

        let poly = Poly::from_coefficients_vec(coefs);
        let poly_comm = kzg_comm_scheme.commit(&poly);

        let z = Scalar::rand(rng);
        let proof = kzg_comm_scheme.prove(&poly, z);

        let is_ok = kzg_comm_scheme.verify(&poly_comm, &proof, &z);

        assert!(is_ok)
    }
}
