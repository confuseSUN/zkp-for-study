use std::marker::PhantomData;

use ark_ff::Field;

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

    use ark_ff::{One, Zero};
    use sample_field::F101;

    #[test]
    pub fn test_boolean_hypercube() {
        let r = BooleanHypercube::new(3)
            .into_iter()
            .collect::<Vec<Vec<F101>>>();

        let one = F101::one();
        let zero = F101::zero();

        let expect = vec![
            vec![zero, zero, zero],
            vec![zero, zero, one],
            vec![zero, one, zero],
            vec![zero, one, one],
            vec![one, zero, zero],
            vec![one, zero, one],
            vec![one, one, zero],
            vec![one, one, one],
        ];

        assert_eq!(r, expect)
    }
}
