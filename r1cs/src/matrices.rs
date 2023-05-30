use ark_ff::Field;

use crate::linear_combination::VarIndex;

#[derive(Clone, Debug, PartialEq)]
pub struct SparseMatrices<F: Field>(pub Vec<Vec<(VarIndex, F)>>);

#[derive(Clone, Debug, PartialEq)]
pub struct DenseMatrices<F: Field>(pub Vec<Vec<F>>);

#[derive(Clone, Debug, PartialEq)]
pub struct R1CSMatrices<M: Matrices> {
    pub a: M,
    pub b: M,
    pub c: M,
}

pub trait Matrices {}

impl<F: Field> Matrices for SparseMatrices<F> {}

impl<F: Field> Matrices for DenseMatrices<F> {}
