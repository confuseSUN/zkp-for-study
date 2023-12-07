use std::time::Instant;

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

#[test]
fn test_comparison() {
    let coeffs = vec![Fr::from(100); 1048576];
    let domain = evaluation_domain(coeffs.len());

    /*   Custom Radix-2 FFT  */
    let start = Instant::now();
    let values = radix2_fft(&coeffs, &domain);
    println!("radix2 fft : {:02?}", start.elapsed());

    /*   Arkwork‘s FFT  */
    let domain = Radix2EvaluationDomain::<Fr>::new(coeffs.len()).unwrap();
    let start = Instant::now();
    let expected = domain.fft(&coeffs);
    println!("arkwork's fft : {:02?}", start.elapsed());

    assert_eq!(values, expected);
}
