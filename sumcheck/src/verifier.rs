use ark_ff::{Field, Zero};
use ark_poly::{univariate::SparsePolynomial, Polynomial};
use ark_std::test_rng;

use crate::SumCheck;

pub struct Verifier<F: Field, S: SumCheck<F> + Default> {
    /// g is the multivariate polynomial to be proved
    g: S,
    // univariate polynolial sent from prover
    intermediate_univariate_poly: SparsePolynomial<F>,
    // random number selected by verifier during sumcheck protocol
    r: Vec<F>,
    /// the sum of g’s evaluations over the Boolean hypercube
    sum: F,
}

impl<F: Field, S: SumCheck<F> + Default> Verifier<F, S> {
    pub fn new(g: S, sum: F) -> Self {
        Self {
            g,
            sum,
            r: vec![],
            intermediate_univariate_poly: SparsePolynomial::zero(),
        }
    }

    pub fn receive_msg(&mut self, j: usize, g_i: SparsePolynomial<F>) -> F {
        let mut rng = test_rng();
        let r = F::rand(&mut rng);

        let eval_0 = g_i.evaluate(&F::zero());
        let eval_1 = g_i.evaluate(&F::one());
        let eval_sum = eval_0.add(eval_1);

        if j == 0 {
            // the first round
            assert_eq!(eval_sum, self.sum);
        } else if 0 < j && j < self.g.num_round() - 1 {
            // the jth round
            let eval_r = self
                .intermediate_univariate_poly
                .evaluate(&self.r.last().unwrap());
            assert_eq!(eval_sum, eval_r);
        } else if j == self.g.num_round() - 1 {
            // the last round
            let eval_r = self
                .intermediate_univariate_poly
                .evaluate(&self.r.last().unwrap());
            assert_eq!(eval_sum, eval_r);

            self.r.push(r.clone());

            let eval_full = self.g.evaluate(&self.r);
            let eval_full_g = g_i.evaluate(&r);
            assert_eq!(eval_full, eval_full_g);
        }

        self.intermediate_univariate_poly = g_i;
        self.r.push(r.clone());

        r
    }
}
