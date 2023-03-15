use ark_ff::{Field, Zero};

use ark_poly::{
    multivariate::{self, SparseTerm, Term},
    univariate, DenseMVPolynomial, Polynomial,
};
use ark_std::{
    fmt::Debug,
    hash::Hash,
    ops::{Add, AddAssign, Neg},
    vec::Vec,
};
use boolean_hypercube::BooleanHypercube;

pub mod boolean_hypercube;

pub trait SumCheck<F: Field>: Clone + Debug + Hash + PartialEq + Eq + Add + Neg + Zero {
    /// Evaluates `self` at the given the vector `point` in slice.
    /// If the number of variables does not match, return `None`.
    fn evaluate(&self, point: &[F]) -> Option<F>;

    /// Reduce the number of variables of `self` by fixing the
    /// `partial_point.len()` variables at `partial_point`.
    fn fix_variables(&self, partial_point: &[F]) -> Self;

    /// Returns a list of evaluations over the domain, which is the boolean
    /// hypercube.
    fn to_evaluations(&self) -> Vec<F>;

    /// Returns univariate polynomial, which is the boolean
    /// hypercube except first variable.
    fn to_univariate(&self) -> univariate::SparsePolynomial<F>;
}

impl<F: Field> SumCheck<F> for multivariate::SparsePolynomial<F, SparseTerm> {
    fn evaluate(&self, point: &[F]) -> Option<F> {
        Some(Polynomial::evaluate(self, &point.to_vec()))
    }

    fn fix_variables(&self, partial_point: &[F]) -> Self {
        let num_vars = self.num_vars();
        let fix_num_vars = partial_point.len();
        let mut res = Self::zero();

        for (coeff, term) in &self.terms {
            let mut new_term_vec = Vec::new();
            let mut eval = *coeff;
            let is_exist_fix_var = term.into_iter().any(|(index, _)| fix_num_vars > *index);

            if is_exist_fix_var {
                for (index, power) in term.into_iter() {
                    if fix_num_vars > *index {
                        eval.mul_assign(partial_point[*index].pow(&[*power as u64]));
                    } else {
                        new_term_vec.push((*index - fix_num_vars, *power))
                    }
                }
            } else {
                new_term_vec = term
                    .into_iter()
                    .map(|(index, power)| (*index - 1, *power))
                    .collect();
            }

            let new_term = SparseTerm::new(new_term_vec);

            let poly = multivariate::SparsePolynomial::<F, SparseTerm> {
                num_vars: num_vars - partial_point.len(),
                terms: vec![(eval, new_term)],
            };

            res.add_assign(&poly);
        }

        res
    }

    fn to_evaluations(&self) -> Vec<F> {
        BooleanHypercube::new(self.num_vars() as u32)
            .into_iter()
            .map(|points| Polynomial::evaluate(self, &points))
            .collect()
    }

    fn to_univariate(&self) -> univariate::SparsePolynomial<F> {
        let mut res = univariate::SparsePolynomial::zero();

        for points in BooleanHypercube::<F>::new((self.num_vars() - 1) as u32).into_iter() {
            let mut full_point = vec![F::one()];
            full_point.extend_from_slice(&points);

            for (coeff, term) in self.terms() {
                let mut eval = term.evaluate(&full_point);
                eval.mul_assign(coeff);

                let power = term
                    .iter()
                    .find(|(index, _)| *index == 0)
                    .map(|(_, power)| *power)
                    .unwrap_or(0);

                let poly = univariate::SparsePolynomial::from_coefficients_slice(&[(power, eval)]);

                res.add_assign(&poly);
            }
        }

        res
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::fields::{Fp64, MontBackend, MontConfig};

    use ark_poly::{
        multivariate::{self, Term},
        univariate::SparsePolynomial,
        DenseMVPolynomial,
    };

    use ark_ff::One;
    use ark_ff::PrimeField;

    use crate::SumCheck;

    #[derive(MontConfig)]
    #[modulus = "101"]
    #[generator = "2"]
    struct FrConfig;

    type Fp101 = Fp64<MontBackend<FrConfig, 1>>;

    #[test]
    fn test_fix_variables() {
        // g(x) = 2 *x_1^3 + x_1 * x_3 + x_2 * x_3
        let g = multivariate::SparsePolynomial::from_coefficients_slice(
            3,
            &[
                (
                    Fp101::from_bigint(2u32.into()).unwrap(),
                    multivariate::SparseTerm::new(vec![(0, 3)]),
                ),
                (
                    Fp101::from_bigint(1u32.into()).unwrap(),
                    multivariate::SparseTerm::new(vec![(0, 1), (2, 1)]),
                ),
                (
                    Fp101::from_bigint(1u32.into()).unwrap(),
                    multivariate::SparseTerm::new(vec![(1, 1), (2, 1)]),
                ),
            ],
        );

        let r1 = Fp101::from(89);
        let g1 = g.fix_variables(&[r1]);
        let expect_g1 = multivariate::SparsePolynomial::from_coefficients_slice(
            2,
            &[
                (
                    Fp101::from_bigint(79u32.into()).unwrap(),
                    multivariate::SparseTerm::new(vec![]),
                ),
                (
                    Fp101::from_bigint(89u32.into()).unwrap(),
                    multivariate::SparseTerm::new(vec![(1, 1)]),
                ),
                (
                    Fp101::from_bigint(1u32.into()).unwrap(),
                    multivariate::SparseTerm::new(vec![(0, 1), (1, 1)]),
                ),
            ],
        );
        assert_eq!(g1, expect_g1);

        let r2 = Fp101::from(1);
        let g2 = g1.fix_variables(&[r2]);
        let expect_g2 = multivariate::SparsePolynomial::from_coefficients_slice(
            1,
            &[
                (
                    Fp101::from_bigint(79u32.into()).unwrap(),
                    multivariate::SparseTerm::new(vec![]),
                ),
                (
                    Fp101::from_bigint(90u32.into()).unwrap(),
                    multivariate::SparseTerm::new(vec![(0, 1)]),
                ),
            ],
        );
        assert_eq!(g2, expect_g2);
    }

    #[test]
    pub fn test_evalutions() {
        // g(x) = 2 *x_1^3 + x_1 * x_3 + x_2 * x_3
        let g = multivariate::SparsePolynomial::from_coefficients_slice(
            3,
            &[
                (
                    Fp101::from_bigint(2u32.into()).unwrap(),
                    multivariate::SparseTerm::new(vec![(0, 3)]),
                ),
                (
                    Fp101::from_bigint(1u32.into()).unwrap(),
                    multivariate::SparseTerm::new(vec![(0, 1), (2, 1)]),
                ),
                (
                    Fp101::from_bigint(1u32.into()).unwrap(),
                    multivariate::SparseTerm::new(vec![(1, 1), (2, 1)]),
                ),
            ],
        );

        let evals = g.to_evaluations();
        let sum: Fp101 = evals.iter().sum();
        assert_eq!(sum, Fp101::from(12))
    }

    #[test]
    pub fn test_to_univariate() {
        // g(x) = 2 *x_1^3 + x_1 * x_3 + x_2 * x_3
        let g = multivariate::SparsePolynomial::from_coefficients_slice(
            3,
            &[
                (
                    Fp101::from_bigint(2u32.into()).unwrap(),
                    multivariate::SparseTerm::new(vec![(0, 3)]),
                ),
                (
                    Fp101::from_bigint(1u32.into()).unwrap(),
                    multivariate::SparseTerm::new(vec![(0, 1), (2, 1)]),
                ),
                (
                    Fp101::from_bigint(1u32.into()).unwrap(),
                    multivariate::SparseTerm::new(vec![(1, 1), (2, 1)]),
                ),
            ],
        );

        let uni_poly = g.to_univariate();

        let coeffs = vec![(0, Fp101::one()), (1, Fp101::from(2)), (3, Fp101::from(8))];
        let expect = SparsePolynomial::from_coefficients_vec(coeffs);

        assert_eq!(uni_poly, expect)
    }
}
