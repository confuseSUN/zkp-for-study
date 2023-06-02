use ark_ff::Field;

use crate::linear_combination::LinearCombination;
use crate::linear_combination::Sign;
use crate::linear_combination::Variable;
use crate::linear_combination::Variable::*;
use crate::matrices::R1CSMatrices;
use crate::matrices::SparseMatrices;
use crate::R1CS;

pub struct ConstraintSystem<F: Field> {
    // The number of constraint.
    pub num_constraint: usize,
    // The number of instance.
    pub num_instance: usize,
    // The number of witness.
    pub num_witness: usize,
    // The witnes.
    pub witness: Vec<F>,
    // The instance.
    pub instance: Vec<F>,
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
            num_instance: 1,
            num_witness: 0,
            num_constraint: 0,
            witness: vec![],
            instance: vec![F::one()],
            a: vec![],
            b: vec![],
            c: vec![],
        }
    }

    pub fn get_and_clear_instance_witness(&mut self) -> Vec<F> {
        let mut instance_witness = self.instance.clone();
        instance_witness.splice(1..1, self.witness.clone());
        self.witness.clear();
        instance_witness
    }

    pub fn is_satisfied(&self) -> bool {
        for ((a, b), c) in self.a.iter().zip(self.b.iter()).zip(self.c.iter()) {
            let left = a.eval(&self.witness, &self.instance);
            let right = b.eval(&self.witness, &self.instance);
            let out = c.eval(&self.witness, &self.instance);

            if left * right != out {
                return false;
            }
        }

        true
    }

    pub fn get_num_instance_and_witness(&self) -> usize {
        self.num_witness + self.num_instance
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
                    Witness(i) => (F::one(), *i + 1),
                    Instance(i) => (F::one(), *i + self.num_witness),
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
                    Witness(i) => (F::one(), *i + 1),
                    Instance(i) => (F::one(), *i + self.num_witness),
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
                    Witness(i) => (F::one(), *i + 1),
                    Instance(i) => (F::one(), *i + self.num_witness),
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
            num_instance_witness: self.get_num_instance_and_witness(),
        }
    }
}

impl<F: Field> R1CS<F> for ConstraintSystem<F> {
    fn new_witness(&mut self, value: F) -> Variable<F> {
        let index = self.num_witness;
        self.num_witness += 1;
        self.witness.push(value);

        Witness(index)
    }

    fn new_instance(&mut self, value: F) -> Variable<F> {
        let index = self.num_instance;
        self.num_instance += 1;
        self.instance.push(value);

        Instance(index)
    }

    fn constrain(
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
}

#[cfg(test)]
mod tests {
    use crate::{
        matrices::{DenseMatrices, R1CSMatrices},
        Circuit, R1CS,
    };
    use ark_ff::Zero;
    use ark_ff::{Field, One};
    use sample_field::BN254Fr;

    use super::ConstraintSystem;

    /// The CubicCircuit, which derived from [Vitalik Buterin](https://vitalik.ca/general/2016/12/10/qap.html),
    /// defines the equation y = x^3 + x + 5.
    struct CubicCircuit<F> {
        pub input: F,
    }

    impl<F: Field> CubicCircuit<F> {
        fn new(input: F) -> Self {
            CubicCircuit { input }
        }
    }

    impl<F: Field> Circuit<F> for CubicCircuit<F> {
        fn synthesize<R: R1CS<F>>(&self, cs: &mut R) {
            let one = F::one();
            let five = F::from(5u8);

            let x = cs.new_witness(self.input);
            let sym1 = cs.new_witness(self.input * self.input);
            let y = cs.new_witness(self.input * self.input * self.input);
            let sym2 = cs.new_witness(self.input * self.input * self.input + self.input);
            let out = cs.new_instance(self.input * self.input * self.input + self.input + five);

            cs.constrain(x.into(), x.into(), sym1.into());
            cs.constrain(sym1.into(), x.into(), y.into());
            cs.constrain(y + x, one.into(), sym2.into());
            cs.constrain(sym2 + five, one.into(), out.into());
        }
    }

    #[test]
    fn test_cubic_circuit() {
        let mut cs = ConstraintSystem::new();
        let circuit = CubicCircuit::new(BN254Fr::from(3));
        circuit.synthesize(&mut cs);

        assert!(cs.is_satisfied());

        let witness = cs.get_and_clear_instance_witness();
        let sparse_martices = cs.to_sparse_matrices();
        let dense_matrices: R1CSMatrices<DenseMatrices<BN254Fr>> = sparse_martices.clone().into();

        let (a, b, c) = get_test_matrices();

        assert_eq!(dense_matrices.a, a);
        assert_eq!(dense_matrices.b, b);
        assert_eq!(dense_matrices.c, c);

        assert!(sparse_martices.verify(&witness));
        assert!(dense_matrices.verify(&witness));
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
            [zero, zero, one, zero, zero, zero].to_vec(),
            [zero, one, zero, one, zero, zero].to_vec(),
            [five, zero, zero, zero, one, zero].to_vec(),
        ];
        let b = vec![
            [zero, one, zero, zero, zero, zero].to_vec(),
            [zero, one, zero, zero, zero, zero].to_vec(),
            [one, zero, zero, zero, zero, zero].to_vec(),
            [one, zero, zero, zero, zero, zero].to_vec(),
        ];
        let c = vec![
            [zero, zero, one, zero, zero, zero].to_vec(),
            [zero, zero, zero, one, zero, zero].to_vec(),
            [zero, zero, zero, zero, one, zero].to_vec(),
            [zero, zero, zero, zero, zero, one].to_vec(),
        ];

        (DenseMatrices(a), DenseMatrices(b), DenseMatrices(c))
    }
}
