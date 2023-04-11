use ark_ff::Field;
use constants::{CONSTANTS, MDS};
use num_bigint::BigUint;

pub mod constants;
pub enum RoundType {
    PartialRound,
    FullRound,
}

pub struct Poseidon<F: Field> {
    // The number of partial rounds
    rp: Vec<usize>,
    // The number of full rounds
    rf: usize,
    // The MDS matrix
    mds: Vec<Vec<Vec<F>>>,
    // The round constants
    constants: Vec<Vec<F>>,
    // The input data
    state: Vec<F>,
    // The index of current the MDS matrix,the partial rounds and the round constants
    index: usize,
}

impl<F: Field + From<BigUint>> Poseidon<F> {
    pub fn new() -> Self {
        let mds = MDS
            .iter()
            .map(|x| {
                x.iter()
                    .map(|y| {
                        y.iter()
                            .map(|z| {
                                let k: BigUint = z.parse().unwrap();
                                F::from(k)
                            })
                            .collect::<Vec<F>>()
                    })
                    .collect::<Vec<Vec<F>>>()
            })
            .collect::<Vec<Vec<Vec<F>>>>();

        let constants = CONSTANTS
            .iter()
            .map(|x| {
                x.iter()
                    .map(|y| {
                        let k: BigUint = y.parse().unwrap();
                        F::from(k)
                    })
                    .collect::<Vec<F>>()
            })
            .collect::<Vec<Vec<F>>>();

        Self {
            rp: vec![56, 57, 56, 60, 60, 63, 64, 63, 60, 66, 60, 65],
            rf: 8,
            mds,
            constants,
            state: vec![],
            index: 0,
        }
    }

    pub fn sbox(&mut self, ith_round: usize) {
        let round_type = self.round_type(ith_round);
        match round_type {
            RoundType::PartialRound => {
                let mut sq = self.state.get(0).unwrap().square();
                sq.square_in_place();
                self.state.get_mut(0).unwrap().mul_assign(&sq)
            }
            RoundType::FullRound => self.state.iter_mut().for_each(|x| {
                let mut sq = x.square();
                sq.square_in_place();
                x.mul_assign(&sq)
            }),
        }
    }

    pub fn arc(&mut self, ith_round: usize) {
        let constant = self.constants.get(self.index).unwrap();
        let constant = &constant
            [ith_round * self.state.len()..ith_round * self.state.len() + self.state.len()];
        self.state
            .iter_mut()
            .zip(constant.iter())
            .for_each(|(s, c)| s.add_assign(c))
    }

    pub fn mix(&mut self) {
        let mds = self.mds.get(self.index).unwrap();
        let mut new_state: Vec<F> = vec![F::zero(); self.state.len()];
        for i in 0..self.state.len() {
            for j in 0..self.state.len() {
                new_state[i].add_assign(&self.state[j].mul(mds[i][j]));
            }
        }

        self.state = new_state
    }

    pub fn permutation(&mut self, input_data: &[F]) -> Result<F, String> {
        if input_data.len() < 1
            || input_data.len() > self.mds.last().unwrap().last().unwrap().len() - 1
        {
            return Err(format!(
                "the length of input data cannot less than 1 and great than {}",
                self.mds.last().unwrap().last().unwrap().len() - 1
            ));
        }

        let mut state = vec![F::zero(); input_data.len() + 1];
        state[1..].copy_from_slice(input_data);

        self.state = state;
        self.index = self.state.len() - 2;

        for i in 0..(self.rf + self.rp[self.index]) {
            self.arc(i);
            self.sbox(i);
            self.mix();
        }

        Ok(self.state[0])
    }

    pub fn round_type(&self, ith_round: usize) -> RoundType {
        if self.rf / 2 <= ith_round && ith_round < self.rf / 2 + self.rp[self.index] {
            RoundType::PartialRound
        } else {
            RoundType::FullRound
        }
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::One;
    use sample_field::BN254Fr;

    use crate::Poseidon;

    #[test]
    // Test vector https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/test_vectors.txt
    fn poseidon_2() {
        let one = BN254Fr::one();
        let two = one + one;

        let input = vec![one, two];

        let mut poseidon = Poseidon::new();

        let r = poseidon.permutation(&input).unwrap();

        assert_eq!(
            r.to_string(),
            "7853200120776062878684798364095072458815029376092732009249414926327459813530"
                .to_string()
        )
    }

    #[test]
    // Test vector https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/test_vectors.txt
    fn poseidon_4() {
        let one = BN254Fr::one();
        let two = one + one;
        let three = one + two;
        let four = two + two;

        let input = vec![one, two, three, four];

        let mut poseidon = Poseidon::new();

        let r = poseidon.permutation(&input).unwrap();

        assert_eq!(
            r.to_string(),
            "18821383157269793795438455681495246036402687001665670618754263018637548127333"
                .to_string()
        )
    }
}
