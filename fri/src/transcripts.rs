use ark_ff::PrimeField;
use ark_std::rand::SeedableRng;
use merlin::Transcript;
use rand_chacha::ChaChaRng;

pub trait GlobalTranscript<F: PrimeField> {
    fn append_scalar(&mut self, scalar: &F);

    fn append_scalars(&mut self, scalar: &[F]);

    fn get_challenge(&mut self, label: &'static [u8]) -> F;
}

impl<F: PrimeField> GlobalTranscript<F> for Transcript {
    fn get_challenge(&mut self, label: &'static [u8]) -> F {
        let mut buf = [0u8; 32];
        self.challenge_bytes(label, &mut buf);
        F::rand(&mut ChaChaRng::from_seed(buf))
    }

    fn append_scalars(&mut self, scalars: &[F]) {
        for scalar in scalars {
            let mut buf = Vec::new();
            scalar.serialize_uncompressed(&mut buf).unwrap();
            self.append_message(b"append scalar", &buf)
        }
    }

    fn append_scalar(&mut self, scalar: &F) {
        let mut buf = Vec::new();
        scalar.serialize_uncompressed(&mut buf).unwrap();
        self.append_message(b"append scalar", &buf)
    }
}
