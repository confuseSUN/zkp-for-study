use ark_ec::AffineCurve;
use ark_ff::PrimeField;
use ark_std::rand::SeedableRng;
use merlin::Transcript;
use rand_chacha::ChaChaRng;

pub trait GlobalTranscript {
    fn append_scalar<T: PrimeField>(&mut self, scalar: &T);

    fn append_scalars<T: PrimeField>(&mut self, scalar: &[T]);

    fn append_points<G: AffineCurve>(&mut self, point: &G);

    fn get_challenge<T: PrimeField>(&mut self, label: &'static [u8]) -> T;
}

impl GlobalTranscript for Transcript {
    fn get_challenge<T: PrimeField>(&mut self, label: &'static [u8]) -> T {
        let mut buf = [0u8; 32];
        self.challenge_bytes(label, &mut buf);
        T::rand(&mut ChaChaRng::from_seed(buf))
    }

    fn append_scalars<T: PrimeField>(&mut self, scalars: &[T]) {
        for scalar in scalars {
            let mut buf = Vec::new();
            scalar.serialize_unchecked(&mut buf).unwrap();
            self.append_message(b"append scalar", &buf)
        }
    }

    fn append_scalar<T: PrimeField>(&mut self, scalar: &T) {
        let mut buf = Vec::new();
        scalar.serialize_unchecked(&mut buf).unwrap();
        self.append_message(b"append scalar", &buf)
    }

    fn append_points<G: AffineCurve>(&mut self, point: &G) {
        let mut buf = Vec::new();
        point.serialize_unchecked(&mut buf).unwrap();
        self.append_message(b"append point", &buf)
    }
}
