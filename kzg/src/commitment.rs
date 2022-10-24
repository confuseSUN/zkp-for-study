use ark_bls12_381::g1::Parameters as G1_Parameters;
use ark_bls12_381::Bls12_381;
use ark_ec::{
    msm::VariableBaseMSM, short_weierstrass_jacobian::GroupAffine, PairingEngine, ProjectiveCurve,
};
use ark_ff::{One, PrimeField, UniformRand};
use ark_poly::{univariate::DenseOrSparsePolynomial, Polynomial, UVPolynomial};
use ark_std::rand::Rng;
use std::ops::{AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::{srs::SRS, Scalar, G1};

type Poly = ark_poly::polynomial::univariate::DensePolynomial<Scalar>;

pub struct KZGCommitmentScheme<'a>(&'a SRS);

#[derive(Debug, Clone)]
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

    pub fn batch_verify<R: Rng>(
        &self,
        poly_comms: &[G1],
        proofs: &[KZGCommitmentProof],
        zs: &[Scalar],
        rng: &mut R,
    ) -> bool {
        let g1_base = self.0.g1[0];
        let g2_base = self.0.g2[0];
        let g2_r = self.0.g2[1];

        let mut r = Scalar::rand(rng);
        let mut left_0 = G1::default();
        left_0.add_assign(&poly_comms[0]);
        left_0.sub_assign(g1_base.mul(&proofs[0].opening_value.into_repr()));
        left_0.add_assign(proofs[0].comm_h.mul(&zs[0].into_repr()));

        let mut right_0 = G1::default();
        right_0.add_assign(&proofs[0].comm_h);

        for ((poly_comm, proof), z) in poly_comms.iter().zip(proofs.iter()).zip(zs.iter()).skip(1) {
            left_0.add_assign(poly_comm.mul(&r.into_repr()));
            left_0.sub_assign(g1_base.mul(&proof.opening_value.into_repr()).mul(&r.into_repr()));
            left_0.add_assign(proof.comm_h.mul(z.mul(&r).into_repr()));

            right_0.add_assign(&proof.comm_h.mul(&r.into_repr()));
           
            r.mul_assign(&r.clone());
        }

        let left = Bls12_381::pairing(left_0, g2_base);

        let right = Bls12_381::pairing(right_0, g2_r);

        left == right
    }
}

#[cfg(test)]
mod test_kzg {
    use super::{KZGCommitmentProof, KZGCommitmentScheme, Poly, SRS};
    use crate::{Scalar, G1};
    use ark_ff::UniformRand;
    use ark_poly::UVPolynomial;
    use ark_std::{rand::Rng, test_rng};

    #[test]
    fn test_kzg_comm() {
        let max_degree = 20;
        let mut rng = test_rng();
        let mut srs = SRS::new(max_degree, &mut rng);
        let srs_clone = srs.clone();
        let kzg_comm_scheme = KZGCommitmentScheme(&srs_clone);
        let (poly_comm, proof, z) = kzg_comm(&kzg_comm_scheme, max_degree, &mut rng);
        let is_ok = kzg_comm_scheme.verify(&poly_comm, &proof, &z);
        assert!(is_ok);

        srs.update(&mut rng);
        let (poly_comm, proof, z) = kzg_comm(&kzg_comm_scheme, max_degree, &mut rng);
        let is_ok = kzg_comm_scheme.verify(&poly_comm, &proof, &z);
        assert!(is_ok);
    }

    #[ignore = "git lfs disable"]
    #[test]
    fn test_kzg_comm_with_public_srs() {
        let max_degree = 20;
        let mut rng = test_rng();
        let mut pub_srs = SRS::load_from_public_setup_parameters(max_degree);
        let pub_srs_clone = pub_srs.clone();
        let kzg_comm_scheme = KZGCommitmentScheme(&pub_srs_clone);
        let (poly_comm, proof, z) = kzg_comm(&kzg_comm_scheme, max_degree, &mut rng);
        let is_ok = kzg_comm_scheme.verify(&poly_comm, &proof, &z);
        assert!(is_ok);

        pub_srs.update(&mut rng);
        let (poly_comm, proof, z) = kzg_comm(&kzg_comm_scheme, max_degree, &mut rng);
        let is_ok = kzg_comm_scheme.verify(&poly_comm, &proof, &z);
        assert!(is_ok);
    }

    #[test]
    fn test_batch_kzg_comm() {
        let max_degree = 20;
        let mut rng = test_rng();
        let srs = SRS::new(max_degree, &mut rng);
        let kzg_comm_scheme = KZGCommitmentScheme(&srs);
        let (poly_comm1, proof1, z1) = kzg_comm(&kzg_comm_scheme, max_degree, &mut rng);
        let is_ok = kzg_comm_scheme.verify(&poly_comm1, &proof1, &z1);
        assert!(is_ok);

        let (poly_comm2, proof2, z2) = kzg_comm(&kzg_comm_scheme, max_degree, &mut rng);
        let is_ok = kzg_comm_scheme.verify(&poly_comm2, &proof2, &z2);
        assert!(is_ok);

        let is_ok = kzg_comm_scheme.batch_verify(
            &[poly_comm1, poly_comm2],
            &[proof1.clone(), proof2.clone()],
            &[z1, z2],
            &mut rng,
        );
        assert!(is_ok);

        let (poly_comm3, proof3, z3) = kzg_comm(&kzg_comm_scheme, max_degree, &mut rng);
        let is_ok = kzg_comm_scheme.verify(&poly_comm3, &proof3, &z3);
        assert!(is_ok);


        let is_ok = kzg_comm_scheme.batch_verify(
            &[poly_comm1, poly_comm2,poly_comm3],
            &[proof1, proof2, proof3],
            &[z1, z2,z3],
            &mut rng,
        );
        assert!(is_ok);
    }

    fn kzg_comm<R: Rng>(
        kzg_comm_scheme: &KZGCommitmentScheme,
        max_degree: usize,
        rng: &mut R,
    ) -> (G1, KZGCommitmentProof, Scalar) {
        let mut coefs = Vec::new();
        for _ in 0..max_degree + 1 {
            let coef = Scalar::rand(rng);
            coefs.push(coef);
        }

        let poly = Poly::from_coefficients_vec(coefs);
        let poly_comm = kzg_comm_scheme.commit(&poly);

        let z = Scalar::rand(rng);
        let proof = kzg_comm_scheme.prove(&poly, z);

        (poly_comm, proof, z)
    }
}
