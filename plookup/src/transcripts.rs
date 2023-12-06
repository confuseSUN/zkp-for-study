use ark_ec::CurveGroup;
use ark_ff::UniformRand;
use ark_std::rand::SeedableRng;
use merlin::Transcript;
use rand_chacha::ChaChaRng;

pub trait GlobalTranscript<G: CurveGroup> {
    fn append_commitent(&mut self, comm: &G);

    fn get_challenge(&mut self, label: &'static [u8]) -> G::ScalarField;
}

impl<G: CurveGroup> GlobalTranscript<G> for Transcript {
    fn append_commitent(&mut self, comm: &G) {
        let mut buf = Vec::new();
        comm.serialize_uncompressed(&mut buf).unwrap();
        self.append_message(b"append commitment", &buf)
    }

    fn get_challenge(&mut self, label: &'static [u8]) -> G::ScalarField {
        let mut buf = [0u8; 32];
        self.challenge_bytes(label, &mut buf);
        G::ScalarField::rand(&mut ChaChaRng::from_seed(buf))
    }
}
