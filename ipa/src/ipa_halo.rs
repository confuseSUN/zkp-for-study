use ark_ec::{msm::VariableBaseMSM, AffineCurve, ProjectiveCurve};
use ark_ff::PrimeField;
use ark_std::test_rng;
use ark_std::UniformRand;
use merlin::Transcript;

use crate::transcripts;
use crate::transcripts::GlobalTranscript;

pub struct Prover<G: ProjectiveCurve> {
    d: usize,
    G: Vec<G>,
    H: G,
}

impl<G: ProjectiveCurve> Prover<G> {
    pub fn new(d: usize) -> Result<Self, String> {
        if !d.is_power_of_two() {
            return Err("d must be power of two".to_string());
        }

        let mut rng = test_rng();
        let H = G::rand(&mut rng);
        let G = (0..d).into_iter().map(|_| G::rand(&mut rng)).collect();
        Ok(Self { d, G, H })
    }

    pub fn commit(&self, a: &[G::ScalarField]) -> (G::ScalarField, G) {
        let mut rng = test_rng();
        let blind = G::ScalarField::rand(&mut rng);
        let mut comm = inner_product_of_scalars_points(a, &self.G);
        comm.add_assign(&self.H.into_affine().mul(blind));
        (blind, comm)
    }

    pub fn commit_with_blind(&self, a: &[G::ScalarField], blind: &G::ScalarField) -> G {
        let mut comm = inner_product_of_scalars_points(a, &self.G);
        comm.add_assign(&self.H.into_affine().mul(*blind));
        comm
    }

    pub fn prove(&self, r: &G::ScalarField, a: &[G::ScalarField], b: &[G::ScalarField], U: &G) {
        if a.len() != self.d || b.len() != self.d {}

        let mut rng = test_rng();
        let k = self.d.trailing_zeros();

        let mut n = self.d;

        let mut transcript = Transcript::new(b"ipa");
        transcript.append_u64(b"size", self.d as u64);
        transcript.append_points(&U.into_affine());

        let mut u: Vec<G::ScalarField> = vec![];

        while n != 1 {
            n = n / 2;

            let a_lo = &a[..n];
            let a_hi = &a[n..];
            let b_lo = &b[..n];
            let b_hi = &b[n..];
            let G_lo = &self.G[..n];
            let G_hi = &self.G[n..];

            let l_j = G::ScalarField::rand(&mut rng);
            let r_j = G::ScalarField::rand(&mut rng);

            let L_j = inner_product_of_scalars_points(a_lo, G_hi)
                .add(self.H.into_affine().mul(l_j))
                .add(
                    U.into_affine()
                        .mul(inner_product_of_scalars_scalars(a_lo, b_hi)),
                );
            let R_j = inner_product_of_scalars_points(a_hi, G_lo)
                .add(self.H.into_affine().mul(r_j))
                .add(
                    U.into_affine()
                        .mul(inner_product_of_scalars_scalars(a_hi, b_lo)),
                );

            transcript.append_points(&L_j.into_affine());
            transcript.append_points(&R_j.into_affine());

            let u_j = transcript.get_challenge(b"u_j");

            u.push(u_j);
        }
    }
}

// We can directly calculate the inner product of type scalars * points using MSM
pub fn inner_product_of_scalars_points<G: ProjectiveCurve>(
    scalars: &[G::ScalarField],
    bases: &[G],
) -> G {
    let scalars = scalars
        .iter()
        .map(|s| s.into_repr())
        .collect::<Vec<<G::ScalarField as PrimeField>::BigInt>>();
    let bases = bases.iter().map(|b| b.into_affine()).collect::<Vec<_>>();
    VariableBaseMSM::multi_scalar_mul(&bases, &scalars)
}

pub fn inner_product_of_scalars_scalars<F: PrimeField>(a: &[F], b: &[F]) -> F {
    a.iter().zip(b.iter()).map(|(x, y)| x.mul(y)).sum()
}

pub fn vec_mul<F: PrimeField>(a: &[F], multiple: &F) -> Vec<F> {
    a.iter().map(|x| x.mul(multiple)).collect()
}

pub fn vec_add<F: PrimeField>(a: &[F], b: &[F]) -> Vec<F> {
    a.iter().zip(b.iter()).map(|(x, y)| x.add(y)).collect()
}

#[cfg(test)]
mod tests {
    use std::ops::{Add, Mul};

    use ark_bls12_381::{Fr, G1Affine, G1Projective};
    use ark_ec::{AffineCurve, ProjectiveCurve};
    use ark_ff::UniformRand;
    use ark_std::test_rng;

    use crate::ipa_halo::inner_product_of_scalars_points;

    use super::{inner_product_of_scalars_scalars, vec_add, vec_mul, Prover};

    #[test]
    fn inner_product_test() {
        let a = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];
        let b = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];
        let G = vec![G1Projective::prime_subgroup_generator(); 4];

        let c = inner_product_of_scalars_scalars(&a, &b);
        assert_eq!(c, Fr::from(30));

        let c = inner_product_of_scalars_points(&a, &G);
        let ten = Fr::from(10);
        assert_eq!(c, G1Affine::prime_subgroup_generator().mul(ten));
    }

    #[test]
    fn add_homomorphic_test() {
        let mut rng = test_rng();

        let a = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];
        let b = vec![Fr::from(11), Fr::from(12), Fr::from(13), Fr::from(14)];

        let prover = Prover::<G1Projective>::new(4).unwrap();

        let (a_blind, a_comm) = prover.commit(&a);
        let (b_blind, b_comm) = prover.commit(&b);

        let a_rand = Fr::rand(&mut rng);
        let b_rand = Fr::rand(&mut rng);
        let a_tmp = vec_mul(&a, &a_rand);
        let b_tmp = vec_mul(&b, &b_rand);

        let c_comm = prover.commit_with_blind(
            &vec_add(&a_tmp, &b_tmp),
            &a_rand.mul(&a_blind).add(b_rand.mul(&b_blind)),
        );

        let expect = a_comm
            .into_affine()
            .mul(a_rand)
            .add(b_comm.into_affine().mul(b_rand));

        assert_eq!(c_comm, expect)
    }
}
