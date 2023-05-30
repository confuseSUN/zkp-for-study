use ark_ff::Field;

use crate::linear_combination::LinearCombination;
use crate::linear_combination::Sign;
use crate::linear_combination::Variable;
use crate::linear_combination::Variable::*;
use crate::matrices::R1CSMatrices;
use crate::matrices::SparseMatrices;

pub struct ConstraintSystem<F: Field> {
    // The number of instance and witness.
    pub num_instance_witness: usize,
    // The number of constraint.
    pub num_constraint: usize,
    // The witness and instance.
    pub instance_witness: Vec<F>,
    // The Constraint of A.
    pub a: Vec<LinearCombination<F>>,
    // The Constraint of B.
    pub b: Vec<LinearCombination<F>>,
    // The Constraint of C.
    pub c: Vec<LinearCombination<F>>,
}

impl<F: Field> ConstraintSystem<F> {
    pub fn new() -> Self {
        ConstraintSystem {
            num_instance_witness: 1,
            num_constraint: 0,
            instance_witness: vec![F::one()],
            a: vec![],
            b: vec![],
            c: vec![],
        }
    }

    pub fn new_witness(&mut self, value: F) -> Variable<F> {
        let index = self.num_instance_witness;
        self.num_instance_witness += 1;
        self.instance_witness.push(value);

        Witness(index)
    }

    pub fn new_instance(&mut self, value: F) -> Variable<F> {
        let index = self.num_instance_witness;
        self.num_instance_witness += 1;
        self.instance_witness.push(value);

        Instance(index)
    }

    pub fn constrain(
        &mut self,
        a: LinearCombination<F>,
        b: LinearCombination<F>,
        c: LinearCombination<F>,
    ) {
        self.num_constraint += 1;
        self.a.push(a);
        self.b.push(b);
        self.c.push(c);
    }

    pub fn to_sparse_matrices(&self) -> R1CSMatrices<SparseMatrices<F>> {
        let mut a = vec![];
        let mut b = vec![];
        let mut c = vec![];

        for lc in self.a.iter() {
            let mut tmp = vec![];
            for (var, sign) in lc.terms.iter() {
                let (mut v, index) = match var {
                    Constant(v) => (*v, 0usize),
                    Instance(i) | Witness(i) => (F::one(), *i),
                };

                if let Sign::Negative = sign {
                    v = -v
                }

                tmp.push((index, v))
            }

            a.push(tmp)
        }

        for lc in self.b.iter() {
            let mut tmp = vec![];
            for (var, sign) in lc.terms.iter() {
                let (mut v, index) = match var {
                    Constant(v) => (*v, 0usize),
                    Instance(i) | Witness(i) => (F::one(), *i),
                };

                if let Sign::Negative = sign {
                    v = -v
                }

                tmp.push((index, v))
            }

            b.push(tmp)
        }

        for lc in self.c.iter() {
            let mut tmp = vec![];
            for (var, sign) in lc.terms.iter() {
                let (mut v, index) = match var {
                    Constant(v) => (*v, 0usize),
                    Instance(i) | Witness(i) => (F::one(), *i),
                };

                if let Sign::Negative = sign {
                    v = -v
                }

                tmp.push((index, v))
            }

            c.push(tmp)
        }

        R1CSMatrices {
            a: SparseMatrices(a),
            b: SparseMatrices(b),
            c: SparseMatrices(c),
            num_instance_witness: self.num_instance_witness,
        }
    }

    pub fn get_and_clear_instance_witness(&mut self) -> Vec<F> {
        let instance_witness = self.instance_witness.clone();
        self.instance_witness.clear();

        instance_witness
    }

    pub fn is_satisfied(&self) -> bool {
        for ((a, b), c) in self.a.iter().zip(self.b.iter()).zip(self.c.iter()) {
            let left = a.eval(&self.instance_witness);
            let right = b.eval(&self.instance_witness);
            let out = c.eval(&self.instance_witness);

            if left * right != out {
                return false;
            }
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        linear_combination::Variable::Constant,
        matrices::{DenseMatrices, R1CSMatrices},
    };
    use ark_ff::One;
    use ark_ff::Zero;
    use sample_field::BN254Fr;

    use super::ConstraintSystem;

    /// The test case is sourced from [Vitalik Buterin](https://vitalik.ca/general/2016/12/10/qap.html).
    /// where y = x ^ 3 + x + 5
    #[test]
    fn test_r1cs_1() {
        let mut cs = ConstraintSystem::new();

        let one = BN254Fr::one();
        let three = BN254Fr::from(3u8);
        let five = BN254Fr::from(5u8);
        let nine = BN254Fr::from(9u8);

        let x = cs.new_witness(three);
        let out = cs.new_instance(nine * three + three + five);
        let sym1 = cs.new_witness(nine);
        let y = cs.new_witness(nine * three);
        let sym2 = cs.new_witness(nine * three + three);

        cs.constrain(x.into(), x.into(), sym1.into());
        cs.constrain(sym1.into(), x.into(), y.into());
        cs.constrain(y + x, one.into(), sym2.into());
        cs.constrain(sym2 + five, one.into(), out.into());

        assert!(cs.is_satisfied());

        let sparse_martices = cs.to_sparse_matrices();

        let dense_matrices: R1CSMatrices<DenseMatrices<BN254Fr>> = sparse_martices.into();

        let (a, b, c) = get_test_matrices();

        assert_eq!(dense_matrices.a, a);
        assert_eq!(dense_matrices.b, b);
        assert_eq!(dense_matrices.c, c);
    }

    fn get_test_matrices() -> (
        DenseMatrices<BN254Fr>,
        DenseMatrices<BN254Fr>,
        DenseMatrices<BN254Fr>,
    ) {
        let zero = BN254Fr::zero();
        let one = BN254Fr::one();
        let five = BN254Fr::from(5);

        let a = vec![
            [zero, one, zero, zero, zero, zero].to_vec(),
            [zero, zero, zero, one, zero, zero].to_vec(),
            [zero, one, zero, zero, one, zero].to_vec(),
            [five, zero, zero, zero, zero, one].to_vec(),
        ];
        let b = vec![
            [zero, one, zero, zero, zero, zero].to_vec(),
            [zero, one, zero, zero, zero, zero].to_vec(),
            [one, zero, zero, zero, zero, zero].to_vec(),
            [one, zero, zero, zero, zero, zero].to_vec(),
        ];
        let c = vec![
            [zero, zero, zero, one, zero, zero].to_vec(),
            [zero, zero, zero, zero, one, zero].to_vec(),
            [zero, zero, zero, zero, zero, one].to_vec(),
            [zero, zero, one, zero, zero, zero].to_vec(),
        ];

        (DenseMatrices(a), DenseMatrices(b), DenseMatrices(c))
    }

    // Z = x^3 - xy -x + 2y -8
    //
    // circuit:
    // o1 = x * y
    // o2 = x * x
    // o3 = o2 * x
    // o4 = 2 * y
    // out = (o3 - o1 - x + o4 -1) * 1
    //
    // witness = [1, x ,y, o1, o2, o3, o4,  out]
    #[test]
    fn test_r1cs_2() {
        let mut cs = ConstraintSystem::new();

        let one = BN254Fr::one();
        let two = BN254Fr::from(2u8);
        let three = BN254Fr::from(3u8);
        let seven = BN254Fr::from(7u8);
        let eight = BN254Fr::from(8u8);

        // if  x =3 and y = 7 then
        // witness = [1, 3, 7, 21, 9, 27, 14, 9]
        let x = cs.new_witness(three); // 3
        let y = cs.new_instance(seven); // 7
        let o1 = cs.new_witness(three * seven); // 21
        let o2 = cs.new_witness(three * three); // 9
        let o3 = cs.new_witness(three * three * three); // 27
        let o4 = cs.new_witness(two * seven); //14
        let out =
            cs.new_instance(three * three * three - three * seven - three + two * seven - eight); //9

        cs.constrain(x.into(), y.into(), o1.into());
        cs.constrain(x.into(), x.into(), o2.into());
        cs.constrain(o2.into(), x.into(), o3.into());
        cs.constrain(two.into(), y.into(), o4.into());
        cs.constrain(o3 - o1 - x + o4 - eight, Constant(one).into(), out.into());

        assert!(cs.is_satisfied());

        let witness = cs.get_and_clear_instance_witness();
        let sparse_martices = cs.to_sparse_matrices();
        let dense_matrices: R1CSMatrices<DenseMatrices<BN254Fr>> = sparse_martices.clone().into();

        assert!(sparse_martices.verify(&witness));
        assert!(dense_matrices.verify(&witness));
    }
}
