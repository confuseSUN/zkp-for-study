use ark_ec::{pairing::Pairing, CurveGroup, VariableBaseMSM};
use ark_ff::{One, UniformRand};
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial},
    DenseUVPolynomial, Polynomial,
};
use ark_std::rand::Rng;
use std::ops::{AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::srs::SRS;

pub struct KZGCommitmentScheme<'a, P: Pairing>(pub &'a SRS<P>);

#[derive(Debug, Clone)]
pub struct KZGCommitmentProof<G: CurveGroup> {
    opening_values: Vec<G::ScalarField>,
    comm_h: G,
}

impl<'a, P: Pairing> KZGCommitmentScheme<'a, P> {
    fn max_degree(&self) -> usize {
        self.0.g1.len() - 1
    }

    pub fn commit(&self, poly: &DensePolynomial<P::ScalarField>) -> P::G1 {
        assert!(self.max_degree() >= poly.degree());
        let bases = P::G1::normalize_batch(&self.0.g1[..poly.coeffs().len()]);
        P::G1::msm(&bases, &poly.coeffs()).unwrap()
    }

    pub fn prove(
        &self,
        poly: &DensePolynomial<P::ScalarField>,
        z: &P::ScalarField,
    ) -> KZGCommitmentProof<P::G1> {
        let opening_value = poly.evaluate(z);

        let dividend = poly.sub(&DensePolynomial::from_coefficients_slice(&[opening_value]));
        let divisor = DensePolynomial::from_coefficients_slice(&[z.neg(), P::ScalarField::one()]);

        let (quotient, remainder) =
            DenseOrSparsePolynomial::divide_with_q_and_r(&(&dividend).into(), &(&divisor).into())
                .unwrap();

        assert!(remainder.is_empty());

        let comm_h = self.commit(&quotient);

        KZGCommitmentProof {
            opening_values: vec![opening_value],
            comm_h,
        }
    }

    pub fn batch_prove(
        &self,
        polys: &[DensePolynomial<P::ScalarField>],
        z: &P::ScalarField,
        _challenge: &P::ScalarField,
    ) -> KZGCommitmentProof<P::G1> {
        //let mut v = Fr::one();
        let mut dividend_sum = DensePolynomial::default();
        let mut opening_values = vec![];
        for poly in polys.iter() {
            let opening_value = poly.evaluate(z);
            opening_values.push(opening_value);
            let dividend = poly.sub(&DensePolynomial::from_coefficients_slice(&[opening_value]));
            //  let dividend = dividend.mul(v);
            dividend_sum.add_assign(&dividend);
            //  v.mul_assign(challenge);
        }

        let divisor = DensePolynomial::from_coefficients_slice(&[z.neg(), P::ScalarField::one()]);

        let (quotient, remainder) = DenseOrSparsePolynomial::divide_with_q_and_r(
            &(&dividend_sum).into(),
            &(&divisor).into(),
        )
        .unwrap();

        assert!(remainder.is_empty());

        let comm_h = self.commit(&quotient);

        KZGCommitmentProof {
            opening_values,
            comm_h,
        }
    }

    pub fn verify(
        &self,
        poly_comm: &P::G1,
        proof: &KZGCommitmentProof<P::G1>,
        z: &P::ScalarField,
    ) -> bool {
        let g1_base = self.0.g1[0];
        let g2_base = self.0.g2[0];
        let g2_r = self.0.g2[1];

        let left = P::pairing(poly_comm.sub(g1_base.mul(proof.opening_values[0])), g2_base);
        let right = P::pairing(proof.comm_h, g2_r.sub(g2_base.mul(z)));

        left == right
    }

    pub fn batch_verify<R: Rng>(
        &self,
        poly_comms: &[Vec<P::G1>],
        proofs: &[&KZGCommitmentProof<P::G1>],
        zs: &[P::ScalarField],
        rng: &mut R,
    ) -> bool {
        let g1_base = self.0.g1[0];
        let g2_base = self.0.g2[0];
        let g2_r = self.0.g2[1];

        let mut r = P::ScalarField::rand(rng);
        let mut left_0 = P::G1::default();
        left_0.add_assign(&poly_comms[0].iter().sum::<P::G1>());
        left_0.sub_assign(g1_base.mul(&proofs[0].opening_values.iter().sum::<P::ScalarField>()));
        left_0.add_assign(proofs[0].comm_h.mul(&zs[0]));

        let mut right_0 = P::G1::default();
        right_0.add_assign(&proofs[0].comm_h);
        for ((poly_comm, proof), z) in poly_comms.iter().zip(proofs.iter()).zip(zs.iter()).skip(1) {
            left_0.add_assign(poly_comm.iter().sum::<P::G1>().mul(&r));
            let opening_values_sum = proof.opening_values.iter().sum::<P::ScalarField>();
            left_0.sub_assign(g1_base.mul(&opening_values_sum).mul(&r));
            left_0.add_assign(proof.comm_h.mul(z.mul(&r)));

            right_0.add_assign(&proof.comm_h.mul(&r));

            r.mul_assign(&r.clone());
        }

        let left = P::pairing(left_0, g2_base);

        let right = P::pairing(right_0, g2_r);

        left == right
    }
}

#[cfg(test)]
mod test_kzg {
    use super::{KZGCommitmentProof, KZGCommitmentScheme, SRS};
    use ark_bls12_381::Bls12_381;
    use ark_ec::pairing::Pairing;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_std::UniformRand;
    use ark_std::{rand::Rng, test_rng};

    #[test]
    fn test_kzg_comm() {
        let max_degree = 20;
        let mut rng = test_rng();
        let mut srs = SRS::new(max_degree, &mut rng);
        let srs_clone = srs.clone();
        let kzg_comm_scheme = KZGCommitmentScheme(&srs_clone);
        let (poly_comm, proof, z) =
            kzg_comm::<Bls12_381, _>(&kzg_comm_scheme, max_degree, &mut rng);
        let is_ok = kzg_comm_scheme.verify(&poly_comm, &proof, &z);
        assert!(is_ok);

        srs.update(&mut rng);
        let (poly_comm, proof, z) = kzg_comm(&kzg_comm_scheme, max_degree, &mut rng);
        let is_ok = kzg_comm_scheme.verify(&poly_comm, &proof, &z);
        assert!(is_ok);
    }

    // #[ignore = "git lfs disable"]
    // #[test]
    // fn test_kzg_comm_with_public_srs() {
    //     let max_degree = 20;
    //     let mut rng = test_rng();
    //     let mut pub_srs = SRS::load_from_public_setup_Prameters(max_degree);
    //     let pub_srs_clone = pub_srs.clone();
    //     let kzg_comm_scheme = KZGCommitmentScheme(&pub_srs_clone);
    //     let (poly_comm, proof, z) = kzg_comm(&kzg_comm_scheme, max_degree, &mut rng);
    //     let is_ok = kzg_comm_scheme.verify(&poly_comm, &proof, &z);
    //     assert!(is_ok);

    //     pub_srs.update(&mut rng);
    //     let (poly_comm, proof, z) = kzg_comm(&kzg_comm_scheme, max_degree, &mut rng);
    //     let is_ok = kzg_comm_scheme.verify(&poly_comm, &proof, &z);
    //     assert!(is_ok);
    // }

    #[test]
    fn test_batch_kzg_comm() {
        let max_degree = 20;
        let batch_size = 5;
        let mut rng = test_rng();
        let srs = SRS::new(max_degree, &mut rng);
        let kzg_comm_scheme = KZGCommitmentScheme(&srs);
        let (poly_comm1, proof1, z1) =
            kzg_comm::<Bls12_381, _>(&kzg_comm_scheme, max_degree, &mut rng);
        let is_ok = kzg_comm_scheme.verify(&poly_comm1, &proof1, &z1);
        assert!(is_ok);

        let (poly_comm2, proof2, z2) = kzg_comm(&kzg_comm_scheme, max_degree, &mut rng);
        let is_ok = kzg_comm_scheme.verify(&poly_comm2, &proof2, &z2);
        assert!(is_ok);

        let is_ok = kzg_comm_scheme.batch_verify(
            &[vec![poly_comm1], vec![poly_comm2]],
            &[&proof1, &proof2],
            &[z1, z2],
            &mut rng,
        );
        assert!(is_ok);

        let (poly_comm3, proof3, z3) = kzg_comm(&kzg_comm_scheme, max_degree, &mut rng);
        let is_ok = kzg_comm_scheme.verify(&poly_comm3, &proof3, &z3);
        assert!(is_ok);

        let is_ok = kzg_comm_scheme.batch_verify(
            &[vec![poly_comm1], vec![poly_comm2], vec![poly_comm3]],
            &[&proof1, &proof2, &proof3.clone()],
            &[z1, z2, z3],
            &mut rng,
        );
        assert!(is_ok);

        let (poly_comm4, proof4, z4) =
            batch_kzg_comm(&kzg_comm_scheme, max_degree, batch_size, &mut rng);
        let is_ok = kzg_comm_scheme.batch_verify(
            &[
                vec![poly_comm1],
                vec![poly_comm2],
                vec![poly_comm3],
                poly_comm4,
            ],
            &[&proof1, &proof2.clone(), &proof3.clone(), &proof4],
            &[z1, z2, z3, z4],
            &mut rng,
        );
        assert!(is_ok);
    }

    fn kzg_comm<P: Pairing, R: Rng>(
        kzg_comm_scheme: &KZGCommitmentScheme<P>,
        max_degree: usize,
        rng: &mut R,
    ) -> (P::G1, KZGCommitmentProof<P::G1>, P::ScalarField) {
        let mut coefs = Vec::new();
        for _ in 0..max_degree + 1 {
            let coef = P::ScalarField::rand(rng);
            coefs.push(coef);
        }

        let poly = DensePolynomial::from_coefficients_vec(coefs);
        let poly_comm = kzg_comm_scheme.commit(&poly);

        let z = P::ScalarField::rand(rng);
        let proof = kzg_comm_scheme.prove(&poly, &z);

        (poly_comm, proof, z)
    }

    fn batch_kzg_comm<P: Pairing, R: Rng>(
        kzg_comm_scheme: &KZGCommitmentScheme<P>,
        max_degree: usize,
        batch_size: usize,
        rng: &mut R,
    ) -> (Vec<P::G1>, KZGCommitmentProof<P::G1>, P::ScalarField) {
        let mut polys = vec![];
        let mut poly_comms = vec![];
        for _ in 0..batch_size {
            let mut coefs = Vec::new();
            for _ in 0..max_degree + 1 {
                let coef = P::ScalarField::rand(rng);
                coefs.push(coef);
            }

            let poly = DensePolynomial::from_coefficients_vec(coefs);
            let poly_comm = kzg_comm_scheme.commit(&poly);

            polys.push(poly);
            poly_comms.push(poly_comm);
        }

        let z = P::ScalarField::rand(rng);
        let challenge = P::ScalarField::rand(rng);
        let proof = kzg_comm_scheme.batch_prove(&polys, &z, &challenge);

        (poly_comms, proof, z)
    }
}
