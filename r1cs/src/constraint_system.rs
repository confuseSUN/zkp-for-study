use ark_ff::Field;

use crate::linear_combination::LinearCombination;
use crate::linear_combination::Sign;
use crate::linear_combination::VarIndex;
use crate::linear_combination::Variable;
use crate::linear_combination::Variable::*;
use crate::matrices::Matrices;
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
        }
    }
}

#[cfg(test)]
mod tests {
    use sample_field::BN254Fr;
    use ark_ff::One;
    use crate::linear_combination::Variable::Constant;

    use super::ConstraintSystem;

    #[test]
    fn test() {
        let mut cs = ConstraintSystem::new();

        let one = BN254Fr::one();
        let three = BN254Fr::from(3u8);
        let five = BN254Fr::from(5u8);
        let nine = BN254Fr::from(9u8);
         
        let x = cs.new_witness(three);
        let out = cs.new_instance(nine * three + three + five);
        let sym1  = cs.new_witness(nine);
        let y = cs.new_witness(nine * three);
        let sym2 = cs.new_witness(nine * three + three);

        cs.constrain(x.into(),x.into(),sym1.into());
        cs.constrain(sym1.into(),x.into(),y.into());
        cs.constrain(y+x,Constant(one).into(),sym2.into());
        cs.constrain(sym2 + Constant(five),Constant(one).into(),out.into());

        let r1cs_martices = cs.to_sparse_matrices();

        println!("{:?}",r1cs_martices.a)
    
    }
}
