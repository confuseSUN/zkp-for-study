use ark_bls12_381::{Fr, G1Projective};
use ark_ff::UniformRand;
use ark_serialize::CanonicalSerialize;
use ark_std::rand::SeedableRng;
use merlin::Transcript;
use rand_chacha::ChaChaRng;

pub trait GlobalTranscript {
    fn append_commitent(&mut self, comm: &G1Projective);

    fn get_challenge(&mut self, label: &'static [u8]) -> Fr;
}

impl GlobalTranscript for Transcript {
    fn append_commitent(&mut self, comm: &G1Projective) {
        let mut buf = Vec::new();
        comm.serialize_uncompressed(&mut buf).unwrap();
        self.append_message(b"append commitment", &buf)
    }

    fn get_challenge(&mut self, label: &'static [u8]) -> Fr {
        let mut buf = [0u8; 32];
        self.challenge_bytes(label, &mut buf);
        Fr::rand(&mut ChaChaRng::from_seed(buf))
    }
}
