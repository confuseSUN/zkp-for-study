use ark_ff::PrimeField;
use rand::{rngs::StdRng, SeedableRng};

pub fn test_rng_helper(seed: [u8; 32]) -> StdRng {
    rand::rngs::StdRng::from_seed(seed)
}

pub fn test_rng_helper_from_scalar<T: PrimeField>(scalar: &T) -> StdRng {
    let mut seed = Vec::new();
    scalar.serialize_compressed(&mut seed).unwrap();

    let seed = seed.try_into().unwrap();
    rand::rngs::StdRng::from_seed(seed)
}

#[cfg(test)]
mod test {
    use super::test_rng_helper;
    use rand::Rng;

    #[test]
    pub fn test() {
        let seed = [
            1, 10, 111, 20, 23, 60, 70, 80, 200, 11, 10, 10, 210, 130, 6, 8, 9, 10, 9, 10, 10, 10,
            10, 0, 0, 0, 0, 0, 0, 0, 0, 255,
        ];

        let mut v1 = vec![];
        let mut rng1 = test_rng_helper(seed.clone());
        for _ in 0..100 {
            let x = rng1.gen_range(0..100);
            v1.push(x);
        }

        let mut v2 = vec![];
        let mut rng2 = test_rng_helper(seed.clone());
        for _ in 0..100 {
            let x = rng2.gen_range(0..100);
            v2.push(x);
        }

        assert!(v1 == v2);
    }
}
