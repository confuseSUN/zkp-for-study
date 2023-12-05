use crate::fri::FRI;
use ark_bls12_381::Fr;
use ark_ff::UniformRand;
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain,
};
use ark_std::test_rng;

#[test]
fn test() {
    // d = p * N -1, p = 1/4
    let degree = 63;
    let expansion_factor = 4;
    let num_colinearity_tests = 4;
    let codeword_length = (degree + 1) * expansion_factor;

    let mut rng = test_rng();
    let mut coefs = Vec::new();
    for _ in 0..degree + 1 {
        let coef = Fr::rand(&mut rng);
        coefs.push(coef);
    }
    let poly = DensePolynomial::from_coefficients_vec(coefs);

    let domain = Radix2EvaluationDomain::<Fr>::new(codeword_length).unwrap();
    let codeword = domain.fft(&poly);
    assert!(codeword.len() == codeword_length);

    let fri = FRI::<Fr>::new(codeword_length, expansion_factor, num_colinearity_tests);

    let proof = fri.prove(&codeword);

    fri.verify(&proof);
}
