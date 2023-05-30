use ark_ff::Field;
use ark_std::ops::Mul;

use crate::linear_combination::VarIndex;

#[derive(Clone, Debug, PartialEq)]
pub struct SparseMatrices<F: Field>(pub Vec<Vec<(VarIndex, F)>>);

impl<'a, F: Field> Mul<&'a [F]> for SparseMatrices<F> {
    type Output = DenseMatrices<F>;

    #[inline]
    fn mul(self, rhs: &[F]) -> Self::Output {
        let mut r: Vec<F> = vec![];

        for x in self.0.iter() {
            let tmp = x.iter().map(|(i, v)| rhs[*i].mul(v)).sum();
            r.push(tmp);
        }

        DenseMatrices(vec![r])
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct DenseMatrices<F: Field>(pub Vec<Vec<F>>);

impl<'a, F: Field> Mul<&'a [F]> for DenseMatrices<F> {
    type Output = DenseMatrices<F>;

    #[inline]
    fn mul(self, rhs: &[F]) -> Self::Output {
        assert!(rhs.len() == self.0[0].len());

        let mut r = Vec::with_capacity(self.0.len());

        for x in self.0.iter() {
            let tmp = x.iter().zip(rhs.iter()).map(|(a, b)| a.mul(b)).sum();
            r.push(tmp);
        }

        DenseMatrices(vec![r])
    }
}

/// element-wise product
impl<F: Field> Mul for DenseMatrices<F> {
    type Output = DenseMatrices<F>;

    #[inline]
    fn mul(self, rhs: Self) -> Self::Output {
        let mut r = vec![];

        for (x, y) in self.0.iter().zip(rhs.0.iter()) {
            let tmp = x.iter().zip(y.iter()).map(|(a, b)| a.mul(b)).collect();
            r.push(tmp)
        }

        DenseMatrices(r)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct R1CSMatrices<M: Matrices> {
    pub a: M,
    pub b: M,
    pub c: M,
    pub num_instance_witness: usize,
}

pub trait Matrices: Clone + for<'a> Mul<&'a [Self::F], Output = DenseMatrices<Self::F>> {
    type F: Field;
}
impl<F: Field> Matrices for SparseMatrices<F> {
    type F = F;
}
impl<F: Field> Matrices for DenseMatrices<F> {
    type F = F;
}

impl<F: Field> From<R1CSMatrices<SparseMatrices<F>>> for R1CSMatrices<DenseMatrices<F>> {
    fn from(v: R1CSMatrices<SparseMatrices<F>>) -> R1CSMatrices<DenseMatrices<F>> {
        let mut a = vec![];
        let mut b = vec![];
        let mut c = vec![];

        for x in v.a.0.iter() {
            let mut tmp = vec![F::zero(); v.num_instance_witness];
            for y in x.iter() {
                tmp[y.0] = y.1
            }
            a.push(tmp);
        }

        for x in v.b.0.iter() {
            let mut tmp = vec![F::zero(); v.num_instance_witness];
            for y in x.iter() {
                tmp[y.0] = y.1
            }
            b.push(tmp);
        }

        for x in v.c.0.iter() {
            let mut tmp = vec![F::zero(); v.num_instance_witness];
            for y in x.iter() {
                tmp[y.0] = y.1
            }
            c.push(tmp);
        }

        R1CSMatrices {
            a: DenseMatrices(a),
            b: DenseMatrices(b),
            c: DenseMatrices(c),
            num_instance_witness: v.num_instance_witness,
        }
    }
}

impl<M: Matrices> R1CSMatrices<M> {
    pub fn verify(&self, witness: &[M::F]) -> bool {
        let left: DenseMatrices<M::F> = self.a.clone().mul(witness);
        let right: DenseMatrices<M::F> = self.b.clone().mul(witness);
        let out: DenseMatrices<M::F> = self.c.clone().mul(witness);

        left * right == out
    }
}
