use ark_ec::{msm::VariableBaseMSM, AffineCurve, ProjectiveCurve};
use ark_ff::batch_inversion;
use ark_ff::Field;
use ark_ff::PrimeField;
use ark_std::test_rng;
use ark_std::UniformRand;
use merlin::Transcript;
use std::ops::Add;
use std::ops::AddAssign;
use std::ops::Mul;

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

    pub fn prove(
        &self,
        a: &[G::ScalarField],
        b: &[G::ScalarField],
        r: &G::ScalarField,
        U: &G,
    ) -> Result<IpAProof<G>, String> {
        if a.len() != self.d || b.len() != self.d {
            return Err("the length of a,b and d must be equal".to_string());
        }

        let mut rng = test_rng();
        let mut n = self.d;

        let mut a_vec = a.to_vec();
        let mut b_vec = b.to_vec();
        let mut G = self.G.to_vec();

        let mut transcript = Transcript::new(b"ipa");
        transcript.append_u64(b"size", self.d as u64);
        transcript.append_points(&U.into_affine());

        let mut u: Vec<G::ScalarField> = vec![];
        let mut L = vec![];
        let mut R = vec![];
        let mut r_blind = r.to_owned();

        while n != 1 {
            n = n / 2;

            let a_lo = &a_vec[..n];
            let a_hi = &a_vec[n..];
            let b_lo = &b_vec[..n];
            let b_hi = &b_vec[n..];
            let G_lo = &G[..n];
            let G_hi = &G[n..];

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

            L.push(L_j);
            R.push(R_j);

            transcript.append_points(&L_j.into_affine());
            transcript.append_points(&R_j.into_affine());
            let u_j: G::ScalarField = transcript.get_challenge(b"u_j");

            u.push(u_j);

            let u_j_inv = u_j.inverse().unwrap();

            a_vec = vec_add(&vec_mul(a_hi, &u_j_inv), &vec_mul(a_lo, &u_j));
            b_vec = vec_add(&vec_mul(b_lo, &u_j_inv), &vec_mul(b_hi, &u_j));
            G = vec_add_point(&vec_mul_point(G_lo, &u_j_inv), &vec_mul_point(G_hi, &u_j));

            let u_j_square = u_j.square();
            let u_j_square_inv = u_j_square.inverse().unwrap();
            r_blind.add_assign(u_j_square.mul(&l_j));
            r_blind.add_assign(u_j_square_inv.mul(&r_j));
        }

        let a_0 = a_vec[0];
        let s = build_s(&u);

        // b_0 = <s,b>
        let b_0 = inner_product_of_scalars_scalars(&s, b);
        assert_eq!(b_0, b_vec[0]);

        // G_0 = <s,G>
        let G_0 = inner_product_of_scalars_points(&s, &self.G);
        assert_eq!(G_0, G[0]);

        // sampling random d, s
        let d = G::ScalarField::rand(&mut rng);
        let s = G::ScalarField::rand(&mut rng);

        // schnorr_R = d * (G + b * U) + s * H
        let schnorr_R = U
            .mul(b_0.into_repr())
            .add(G_0)
            .mul(d.into_repr())
            .add(self.H.mul(s.into_repr()));

        transcript.append_points(&schnorr_R.into_affine());
        let c: G::ScalarField = transcript.get_challenge(b"challenge c");

        let z_1 = a_0.mul(&c).add(&d);
        let z_2 = r_blind.mul(&c).add(&s);

        let schnorr_proof = SchnorrProof::<G> {
            R: schnorr_R,
            z_1,
            z_2,
        };

        Ok(IpAProof {
            L,
            R,
            schnorrProof: schnorr_proof,
            G: self.G.clone(),
            H: self.H.clone(),
        })
    }
}

pub struct IpAProof<G: ProjectiveCurve> {
    L: Vec<G>,
    R: Vec<G>,
    schnorrProof: SchnorrProof<G>,
    G: Vec<G>,
    H: G,
}

pub struct SchnorrProof<G: ProjectiveCurve> {
    R: G,
    z_1: G::ScalarField,
    z_2: G::ScalarField,
}

impl<G: ProjectiveCurve> IpAProof<G> {
    fn verify(&self, d: usize, P: &G, U: &G, b: &[G::ScalarField]) -> bool {
        let mut Q = P.to_owned();

        let mut transcript = Transcript::new(b"ipa");
        transcript.append_u64(b"size", d as u64);
        transcript.append_points(&U.into_affine());
        let u: Vec<G::ScalarField> = self
            .L
            .iter()
            .zip(self.R.iter())
            .map(|(L_j, R_j)| {
                transcript.append_points(&L_j.into_affine());
                transcript.append_points(&R_j.into_affine());
                transcript.get_challenge(b"u_j")
            })
            .collect();

        u.iter()
            .zip(self.L.iter())
            .zip(self.R.iter())
            .for_each(|((u_j, L_j), R_j)| {
                let u_j_square = u_j.square();
                let u_j_square_inv = u_j_square.inverse().unwrap();
                Q.add_assign(L_j.mul(u_j_square.into_repr()));
                Q.add_assign(R_j.mul(u_j_square_inv.into_repr()));
            });

        transcript.append_points(&self.schnorrProof.R.into_affine());
        let c: G::ScalarField = transcript.get_challenge(b"challenge c");

        let s = build_s(&u);
        let b_0 = inner_product_of_scalars_scalars(&s, b);
        let G_0 = inner_product_of_scalars_points(&s, &self.G);

        let left = Q.mul(c.into_repr()).add(&self.schnorrProof.R);
        let right = G_0
            .add(U.mul(b_0.into_repr()))
            .mul(self.schnorrProof.z_1.into_repr())
            .add(self.H.mul(self.schnorrProof.z_2.into_repr()));

        if left == right {
            return true;
        }

        false
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

pub fn vec_mul_point<G: ProjectiveCurve, F: PrimeField>(a: &[G], multiple: &F) -> Vec<G> {
    a.iter().map(|x| x.mul(multiple.into_repr())).collect()
}

pub fn vec_add<F: PrimeField>(a: &[F], b: &[F]) -> Vec<F> {
    a.iter().zip(b.iter()).map(|(x, y)| x.add(y)).collect()
}

pub fn vec_add_point<G: ProjectiveCurve>(a: &[G], b: &[G]) -> Vec<G> {
    a.iter().zip(b.iter()).map(|(x, y)| x.add(y)).collect()
}

// s = (
//   u₁⁻¹ u₂⁻¹ … uₖ⁻¹,
//   u₁   u₂⁻¹ … uₖ⁻¹,
//   u₁⁻¹ u₂   … uₖ⁻¹,
//   u₁   u₂   … uₖ⁻¹,
//   ⋮    ⋮      ⋮
//   u₁   u₂   … uₖ
// )
// Hard-coding is used here, with an input length of 8.
fn build_s<F: PrimeField>(u: &[F]) -> Vec<F> {
    let mut s: Vec<F> = vec![F::one(); 8];

    let mut u_inv = u.to_owned();
    batch_inversion(&mut u_inv);

    s[0] = u_inv[0].mul(&u_inv[1]).mul(&u_inv[2]);
    s[1] = u_inv[0].mul(&u_inv[1]).mul(&u[2]);
    s[2] = u_inv[0].mul(&u[1]).mul(&u_inv[2]);
    s[3] = u_inv[0].mul(&u[1]).mul(&u[2]);
    s[4] = u[0].mul(&u_inv[1]).mul(&u_inv[2]);
    s[5] = u[0].mul(&u_inv[1]).mul(&u[2]);
    s[6] = u[0].mul(&u[1]).mul(&u_inv[2]);
    s[7] = u[0].mul(&u[1]).mul(&u[2]);
    s
}

#[cfg(test)]
mod tests {
    use std::ops::{Add, AddAssign, Mul};

    use ark_bls12_381::{Fr, G1Affine, G1Projective};
    use ark_ec::{AffineCurve, ProjectiveCurve};
    use ark_ff::{One, PrimeField, UniformRand};
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

    #[test]
    fn ipa_test() {
        let d = 8;
        let mut rng = test_rng();

        let mut a = vec![];
        for i in 0..d {
            a.push(Fr::from((i + 1) as u64));
        }

        let prover = Prover::<G1Projective>::new(d).unwrap();
        let (comm_blind, mut comm) = prover.commit(&a);

        // verifier select x and U ;
        let x = Fr::rand(&mut rng);
        let U = G1Projective::rand(&mut rng);
        let mut b = vec![Fr::one()];
        for i in 1..d {
            b.push(b[i - 1].mul(&x))
        }

        let v = inner_product_of_scalars_scalars(&a, &b);
        comm.add_assign(U.mul(v.into_repr()));

        let proof = prover.prove(&a, &b, &comm_blind, &U).unwrap();

        let r = proof.verify(d, &comm, &U, &b);

        assert!(r)
    }
}
