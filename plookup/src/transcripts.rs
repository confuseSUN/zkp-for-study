use ark_ff::UniformRand;
use ark_serialize::CanonicalSerialize;
use ark_std::rand::SeedableRng;
use kzg::{Scalar, G1};
use merlin::Transcript;
use rand_chacha::ChaChaRng;

pub trait GlobalTranscript {
    fn append_commitent(&mut self, comm: &G1);

    fn get_challenge(&mut self, label: &'static [u8]) -> Scalar;
}

impl GlobalTranscript for Transcript {
    fn append_commitent(&mut self, comm: &G1) {
        let mut buf = Vec::new();
        comm.serialize_unchecked(&mut buf).unwrap();
        self.append_message(b"append commitment", &buf)
    }

    fn get_challenge(&mut self, label: &'static [u8]) -> Scalar {
        let mut buf = [0u8; 32];
        self.challenge_bytes(label, &mut buf);
        Scalar::rand(&mut ChaChaRng::from_seed(buf))
    }
}
