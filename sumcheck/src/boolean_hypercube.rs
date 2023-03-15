use std::{fmt::Binary, io::Bytes, marker::PhantomData};

use ark_ff::Field;
use bitvec::slice::BitSlice;

pub struct BooleanHypercube<F: Field> {
    n: u32,
    current: u64,
    _phantom_data: PhantomData<F>,
}

impl<F: Field> BooleanHypercube<F> {
    pub fn new(n: u32) -> Self {
        Self {
            n,
            current: 0,
            _phantom_data: PhantomData,
        }
    }
}

impl<F: Field> Iterator for BooleanHypercube<F> {
    type Item = Vec<F>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current != 2u64.pow(self.n) {
            let mut b = format!("{:b}", self.current);
            while self.n as usize > b.len() {
                b = format!("0{}", b)
            }
            self.current += 1;

            Some(
                b.bytes()
                    .map(|f| match f {
                        49 => F::one(),
                        48 => F::zero(),
                        _ => unreachable!(),
                    })
                    .collect(),
            )
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {
    use super::BooleanHypercube;

    use ark_ff::fields::{Fp64, MontBackend, MontConfig};

    #[derive(MontConfig)]
    #[modulus = "101"]
    #[generator = "2"]
    struct FrConfig;
    type Fp101 = Fp64<MontBackend<FrConfig, 1>>;

    #[test]
    pub fn test_boolean_hypercube() {
        let r = BooleanHypercube::new(3)
            .into_iter()
            .collect::<Vec<Vec<Fp101>>>();

        println!("{:?}", r)
    }
}
