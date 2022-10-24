use ark_bls12_381::g1::Parameters as G1_Parameters;
use ark_bls12_381::Bls12_381;
use ark_ec::{
    msm::VariableBaseMSM, short_weierstrass_jacobian::GroupAffine, PairingEngine, ProjectiveCurve,
};
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

        let g1 = self
            .0
            .g1
            .iter()
            .map(|x| x.into_affine())
            .collect::<Vec<GroupAffine<G1_Parameters>>>();
        let (num_leading_zeros, coeffs) = Self::skip_leading_zeros_and_convert_to_bigints(poly);
        let commitment = VariableBaseMSM::multi_scalar_mul(&g1[num_leading_zeros..], &coeffs);

        commitment
    }

    fn skip_leading_zeros_and_convert_to_bigints<F: PrimeField, P: UVPolynomial<F>>(
        p: &P,
    ) -> (usize, Vec<F::BigInt>) {
        let mut num_leading_zeros = 0;
        while num_leading_zeros < p.coeffs().len() && p.coeffs()[num_leading_zeros].is_zero() {
            num_leading_zeros += 1;
        }
        let coeffs = Self::convert_to_bigints(&p.coeffs()[num_leading_zeros..]);
        (num_leading_zeros, coeffs)
    }

    fn convert_to_bigints<F: PrimeField>(p: &[F]) -> Vec<F::BigInt> {
        let coeffs = ark_std::cfg_iter!(p)
            .map(|s| s.into_repr())
            .collect::<Vec<_>>();
        coeffs
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
