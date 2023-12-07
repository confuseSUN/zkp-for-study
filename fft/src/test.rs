use ark_bls12_381::Fr;
use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain,
};

use crate::{domain::evaluation_domain, fft::radix2_fft};

fn native_eval<F: PrimeField>(coeffs: &[F]) -> Vec<F> {
    let domain = Radix2EvaluationDomain::new(coeffs.len()).unwrap();
    let poly = DensePolynomial::from_coefficients_slice(coeffs);
    let r = poly.evaluate_over_domain(domain);
    r.evals
}

#[test]
fn test_fft() {
    let coeffs = vec![
        Fr::from(1),
        Fr::from(2),
        Fr::from(3),
        Fr::from(4),
        Fr::from(5),
        Fr::from(6),
        Fr::from(7),
        Fr::from(8),
    ];
    let domain = evaluation_domain(coeffs.len());
    let values = radix2_fft(&coeffs, &domain);

    let expected = native_eval(&coeffs);

    assert_eq!(values, expected);
}
