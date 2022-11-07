use std::ops::Mul;

use ark_poly::Radix2EvaluationDomain;
use ark_std::rand::Rng;
use merlin::Transcript;
use study_kzg::{
    commitment::{KZGCommitmentProof, KZGCommitmentScheme},
    Scalar, G1,
};

use crate::transcripts::GlobalTranscript;

pub struct PlookUpProof {
    pub f_comm: G1,
    pub t_comm: G1,
    pub h1_comm: G1,
    pub h2_comm: G1,
    pub r_comm: G1,
    pub z_comm: G1,
    pub opening_witness_zeta: KZGCommitmentProof,
    pub opening_witness_zeta_omega: KZGCommitmentProof,
    pub domain: Radix2EvaluationDomain<Scalar>,
}

impl PlookUpProof {
    pub fn verify<R: Rng>(&self, kzg_comm_scheme: &KZGCommitmentScheme, rng: &mut R) -> bool {
        let mut transcript = Transcript::new(b"plookup");
        transcript.append_u64(b"size", self.domain.size as u64);
        transcript.append_commitent(&self.t_comm);
        transcript.append_commitent(&self.f_comm);
        transcript.append_commitent(&self.h1_comm);
        transcript.append_commitent(&self.h2_comm);
        let _beta = transcript.get_challenge(b"beta");
        let _gamma = transcript.get_challenge(b"gamma");
        transcript.append_commitent(&self.z_comm);

        let zeta = transcript.get_challenge(b"zeta");
        let zeta_omega = zeta.mul(self.domain.group_gen);

        transcript.append_commitent(&self.r_comm);
        let _v = transcript.get_challenge(b"v");

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
