use std::ops::Mul;

use ark_ec::{pairing::Pairing, CurveGroup};
use ark_poly::Radix2EvaluationDomain;
use ark_std::rand::Rng;
use kzg::commitment::{KZGCommitmentProof, KZGCommitmentScheme};
use merlin::Transcript;

use crate::transcripts::GlobalTranscript;

pub struct PlookUpProof<G: CurveGroup> {
    pub f_comm: G,
    pub t_comm: G,
    pub h1_comm: G,
    pub h2_comm: G,
    pub r_comm: G,
    pub z_comm: G,
    pub opening_witness_zeta: KZGCommitmentProof<G>,
    pub opening_witness_zeta_omega: KZGCommitmentProof<G>,
    pub domain: Radix2EvaluationDomain<G::ScalarField>,
}

impl<G: CurveGroup> PlookUpProof<G> {
    pub fn verify<P: Pairing<G1 = G, ScalarField = G::ScalarField>, R: Rng>(
        &self,
        kzg_comm_scheme: &KZGCommitmentScheme<P>,
        rng: &mut R,
    ) -> bool {
        let mut transcript = Transcript::new(b"plookup");
        transcript.append_u64(b"size", self.domain.size as u64);
        transcript.append_commitent(&self.t_comm);
        transcript.append_commitent(&self.f_comm);
        transcript.append_commitent(&self.h1_comm);
        transcript.append_commitent(&self.h2_comm);
        let _beta = <Transcript as GlobalTranscript<G>>::get_challenge(&mut transcript, b"beta");
        let _gamma = <Transcript as GlobalTranscript<G>>::get_challenge(&mut transcript, b"gamma");
        transcript.append_commitent(&self.z_comm);

        let zeta = <Transcript as GlobalTranscript<G>>::get_challenge(&mut transcript, b"zeta");
        let zeta_omega = zeta.mul(&self.domain.group_gen);

        transcript.append_commitent(&self.r_comm);
        let _v = <Transcript as GlobalTranscript<G>>::get_challenge(&mut transcript, b"v");

        let is_ok = kzg_comm_scheme.batch_verify(
            &[
                vec![self.f_comm, self.t_comm, self.h1_comm, self.r_comm],
                vec![self.t_comm, self.h1_comm, self.h2_comm, self.z_comm],
            ],
            &[&self.opening_witness_zeta, &self.opening_witness_zeta_omega],
            &[zeta, zeta_omega],
            rng,
        );

        is_ok
    }
}
