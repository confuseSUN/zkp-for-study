use ark_ff::Field;
use ark_poly::univariate::SparsePolynomial;

use crate::SumCheck;

pub struct Prover<F: Field, S: SumCheck<F> + Default> {
    /// g is the multivariate polynomial to be proved
    g: S,
    // intermediate_g is the intermediate result of g executing `fix_variable`
    intermediate_g: S,
    /// the sum of g’s evaluations over the Boolean hypercube
    sum: F,
}

impl<F: Field, S: SumCheck<F> + Default> Prover<F, S> {
    pub fn new(g: S) -> Self {
        let sum = g.to_evaluations().iter().sum();
        Self {
            g,
            sum,
            intermediate_g: S::default(),
        }
    }

    pub fn get_sum(&self) -> F {
        self.sum
    }

    pub fn start_round(&mut self, j: usize, r_i: F) -> SparsePolynomial<F> {
        if j == 0 {
            self.intermediate_g = self.g.clone()
        }

        if j != 0 {
            self.intermediate_g = self.intermediate_g.fix_variables(&[r_i]);
        }

        self.intermediate_g.to_univariate()
    }
}
