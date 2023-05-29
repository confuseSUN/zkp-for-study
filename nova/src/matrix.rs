use ark_ff::Field;
use ark_relations::r1cs::ConstraintMatrices;
use ark_std::ops::Mul;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix<F: Field>(pub Vec<Vec<F>>);

impl<'a, F: Field> Mul<&'a [F]> for Matrix<F> {
    type Output = Matrix<F>;

    #[inline]
    fn mul(self, rhs: &[F]) -> Self::Output {
        assert!(rhs.len() == self.0[0].len());

        let mut r = Vec::with_capacity(self.0.len());

        for x in self.0.iter() {
            let tmp = x.iter().zip(rhs.iter()).map(|(a, b)| a.mul(b)).sum();
            r.push(tmp);
        }

        Matrix(vec![r])
    }
}

/// element-wise product
impl<F: Field> Mul for Matrix<F> {
    type Output = Matrix<F>;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        let mut r = vec![];

        for (x, y) in self.0.iter().zip(rhs.0.iter()) {
            let tmp = x.iter().zip(y.iter()).map(|(a, b)| a.mul(b)).collect();
            r.push(tmp)
        }

        Matrix(r)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct R1CSMatrices<F: Field> {
    pub a: Matrix<F>,
    pub b: Matrix<F>,
    pub c: Matrix<F>,
}

impl<F: Field> R1CSMatrices<F> {
    pub fn verify(&self, witness: &[F]) -> bool {
        let left = self.a.clone().mul(witness);
        let right = self.b.clone().mul(witness);
        let out = self.c.clone().mul(witness);

        left * right == out
    }
}

pub trait SplitCSMatric<F: Field> {
    fn split(&self) -> R1CSMatrices<F>;
}

impl<F: Field> SplitCSMatric<F> for ConstraintMatrices<F> {
    fn split(&self) -> R1CSMatrices<F> {
        let mut a = vec![];
        let mut b = vec![];
        let mut c = vec![];

        for x in self.a.iter() {
            let mut tmp = vec![F::zero(); self.num_instance_variables + self.num_witness_variables];
            for y in x.iter() {
                tmp[y.1] = y.0
            }
            a.push(tmp);
        }

        for x in self.b.iter() {
            let mut tmp = vec![F::zero(); self.num_instance_variables + self.num_witness_variables];
            for y in x.iter() {
                tmp[y.1] = y.0
            }
            b.push(tmp);
        }

        for x in self.c.iter() {
            let mut tmp = vec![F::zero(); self.num_instance_variables + self.num_witness_variables];
            for y in x.iter() {
                tmp[y.1] = y.0
            }
            c.push(tmp);
        }

        R1CSMatrices {
            a: Matrix(a),
            b: Matrix(b),
            c: Matrix(c),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::matrix::SplitCSMatric;
    use ark_ff::One;
    use ark_ff::Zero;
    use ark_relations::{
        lc,
        r1cs::{ConstraintSystem, Variable},
    };
    use ark_test_curves::bls12_381::Fr;

    use super::Matrix;
    use super::R1CSMatrices;

    /// The test case is sourced from [Vitalik Buterin](https://vitalik.ca/general/2016/12/10/qap.html).
    /// where y = x ^ 3 + x + 5
    #[test]
    fn matrix_generation_example() {
        let cs = ConstraintSystem::<Fr>::new_ref();
        let three = Fr::from(3u8);
        let five = Fr::from(5u8);
        let nine = Fr::from(9u8);

        let x = cs.new_witness_variable(|| Ok(three)).unwrap();
        let out = cs
            .new_input_variable(|| Ok(nine * three + three + five))
            .unwrap();
        let sym_1 = cs.new_witness_variable(|| Ok(nine)).unwrap();
        let y = cs.new_witness_variable(|| Ok(nine * three)).unwrap();
        let sym_2 = cs
            .new_witness_variable(|| Ok(nine * three + three))
            .unwrap();

        cs.enforce_constraint(lc!() + x, lc!() + x, lc!() + sym_1)
            .unwrap();
        cs.enforce_constraint(lc!() + sym_1, lc!() + x, lc!() + y)
            .unwrap();
        cs.enforce_constraint(lc!() + y + x, lc!() + Variable::One, lc!() + sym_2)
            .unwrap();
        cs.enforce_constraint(
            lc!() + sym_2 + (five, Variable::One),
            lc!() + Variable::One,
            lc!() + out,
        )
        .unwrap();

        cs.finalize();
        assert!(cs.is_satisfied().unwrap());
        let matrices = cs.to_matrices().unwrap();

        assert_eq!(matrices.a[0], vec![(Fr::one(), 2)]);
        assert_eq!(matrices.b[0], vec![(Fr::one(), 2)]);
        assert_eq!(matrices.c[0], vec![(Fr::one(), 3)]);

        assert_eq!(matrices.a[1], vec![(Fr::one(), 3)]);
        assert_eq!(matrices.b[1], vec![(Fr::one(), 2)]);
        assert_eq!(matrices.c[1], vec![(Fr::one(), 4)]);

        assert_eq!(matrices.a[2], vec![(Fr::one(), 2), (Fr::one(), 4)]);
        assert_eq!(matrices.b[2], vec![(Fr::one(), 0)]);
        assert_eq!(matrices.c[2], vec![(Fr::one(), 5)]);

        assert_eq!(matrices.a[3], vec![(five, 0), (Fr::one(), 5)]);
        assert_eq!(matrices.b[3], vec![(Fr::one(), 0)]);
        assert_eq!(matrices.c[3], vec![(Fr::one(), 1)]);

        let r1cs_matrices = matrices.split();
        let (a, b, c) = get_test_matrices();

        assert_eq!(a.0, r1cs_matrices.a.0);
        assert_eq!(b.0, r1cs_matrices.b.0);
        assert_eq!(c.0, r1cs_matrices.c.0);
    }

    #[test]
    fn test_matrices() {
        let (a, b, c) = get_test_matrices();

        let witness = vec![
            Fr::one(),
            Fr::from(35),
            Fr::from(3),
            Fr::from(9),
            Fr::from(27),
            Fr::from(30),
        ];

        let r1cs_matrices = R1CSMatrices { a, b, c };
        assert!(r1cs_matrices.verify(&witness))
    }

    fn get_test_matrices() -> (Matrix<Fr>, Matrix<Fr>, Matrix<Fr>) {
        let zero = Fr::zero();
        let one = Fr::one();
        let five = Fr::from(5);

        let a = vec![
            [zero, zero, one, zero, zero, zero].to_vec(),
            [zero, zero, zero, one, zero, zero].to_vec(),
            [zero, zero, one, zero, one, zero].to_vec(),
            [five, zero, zero, zero, zero, one].to_vec(),
        ];
        let b = vec![
            [zero, zero, one, zero, zero, zero].to_vec(),
            [zero, zero, one, zero, zero, zero].to_vec(),
            [one, zero, zero, zero, zero, zero].to_vec(),
            [one, zero, zero, zero, zero, zero].to_vec(),
        ];
        let c = vec![
            [zero, zero, zero, one, zero, zero].to_vec(),
            [zero, zero, zero, zero, one, zero].to_vec(),
            [zero, zero, zero, zero, zero, one].to_vec(),
            [zero, one, zero, zero, zero, zero].to_vec(),
        ];

        (Matrix(a), Matrix(b), Matrix(c))
    }
}
