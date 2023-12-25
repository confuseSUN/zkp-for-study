use ark_ff::{FftField, Field};
use ark_poly::Polynomial;
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain,
};
use fibonacci_constraint::{FibonacciConstraint, Label};
use field::Fr;
use std::ops::*;
use utils::vec_div;

use crate::utils::{deep_quotient, deep_quotient_two_points, mix};

pub mod fibonacci_constraint;
pub mod field;
pub mod utils;

fn main() {
    let public_inputs = vec![Fr::from(24), Fr::from(30), Fr::from(28)];

    // 1. The Execution Trace with padding.

    let data_column_1 = vec![24, 30, 54, 84, 78, 15, 29, 50];
    let data_column_2 = vec![30, 54, 84, 41, 2, 77, 21, 36];
    let data_column_3 = vec![54, 84, 41, 28, 71, 17, 92, 33];

    let control_column_1 = vec![1, 0, 0, 0, 0, 0, 0, 0];
    let control_column_2 = vec![0, 1, 1, 1, 0, 0, 0, 0];
    let control_column_3 = vec![0, 0, 0, 1, 0, 0, 0, 0];

    let trace: Vec<Vec<Fr>> = vec![
        data_column_1.iter().map(|x| Fr::from(*x)).collect(),
        data_column_2.iter().map(|x| Fr::from(*x)).collect(),
        data_column_3.iter().map(|x| Fr::from(*x)).collect(),
        control_column_1.iter().map(|x| Fr::from(*x)).collect(),
        control_column_2.iter().map(|x| Fr::from(*x)).collect(),
        control_column_3.iter().map(|x| Fr::from(*x)).collect(),
    ];

    // 2. Constructing Trace Polynomials.

    // println!("The six trace polynomials:");
    let mut trace_polys = Vec::new();
    let domain = Radix2EvaluationDomain::<Fr>::new(trace[0].len()).unwrap();
    for t in trace.iter() {
        let coeffs = domain.ifft(t);
        // println!("{:?}", coeffs);
        trace_polys.push(coeffs);
    }

    // 3. Evaluating the trace polynomial on the larger domain gives a Reed-Solomon encoded trace block.

    //  println!("\n The reedsolomon codes of trace polynomial:");
    let mut trace_reedsolomon_codes = Vec::new();
    let expanded_domain = Radix2EvaluationDomain::<Fr>::new(trace[0].len() * 4).unwrap();
    for coeffs in trace_polys.iter() {
        let evals = expanded_domain.fft(coeffs);
        // println!("{:?} \n", evals);
        trace_reedsolomon_codes.push(evals);
    }

    // 4. ZK Commitments of the Trace Data.

    //  println!("\n The zk commitments of trace polynomial:");
    let mut trace_zk_commitments = Vec::new();
    let expanded_coset_domain = expanded_domain.get_coset(Fr::GENERATOR).unwrap();
    for coeffs in trace_polys.iter() {
        let evals = expanded_coset_domain.fft(coeffs);
        // println!("{:?} \n", evals);
        trace_zk_commitments.push(evals);
    }

    // 5. Constraint Polynomials.

    let constraints = FibonacciConstraint::construct(&trace, &public_inputs, Label::Trace);
    let constraint_reedsolomon_codes = FibonacciConstraint::construct(
        &trace_reedsolomon_codes,
        &public_inputs,
        Label::ReedsolomonOrZkCommitment,
    );
    let constraint_zk_commitments = FibonacciConstraint::construct(
        &trace_zk_commitments,
        &public_inputs,
        Label::ReedsolomonOrZkCommitment,
    );
    // println!("constraints \n:{:?}", constraints);
    // println!("constraint reedsolomon codes:\n{:?}", constraint_reedsolomon_codes);
    // println!("constraint zk commitments:\n{:?}", constraint_zk_commitments);

    // 6. Mixing Constraint Polynomials.

    let mix_alpha = Fr::from(3);
    let _mixed_constraint_column = constraints.mix(&mix_alpha);
    let _mixed_constraint_reedsolomon_code = constraint_reedsolomon_codes.mix(&mix_alpha);
    let mixed_constraint_zk_commitment = constraint_zk_commitments.mix(&mix_alpha);
    // println!("mixed constraint :{:?} \n", mixed_constraint_zk_commitment);

    // let _constraint_polys = constraint_reedsolomon_codes.ifft(&expanded_domain);
    // let _mixed_constraint_polys = expanded_domain.ifft(&mixed_constraint_reedsolomon_code);
    // println!("constraint polys:\n{:?}",constraint_polys);
    // println!("mixed constraint polys:{:?} \n", mixed_constraint_polys);

    //  7. Constructing the validity polynomial by dividing the mixed constraint polynomial by the zeros polynomial.

    let mut z_x = vec![Fr::ZERO; 9]; // domain.vanishing_polynomial();
    z_x[0] = Fr::ONE.neg();
    z_x[8] = Fr::ONE;
    let _z_x_reedsolomon_code = expanded_domain.fft(&z_x);
    let z_x_zk_commitment = expanded_coset_domain.fft(&z_x);
    // println!("reedsolomon code of z_x:{:?}", z_x_reedsolomon_code);
    // println!("zk commitment of z_x:{:?}", z_x_zk_commitment);

    let validity_poly_evals = vec_div(&mixed_constraint_zk_commitment, &z_x_zk_commitment);
    let validity_poly_coeffs = expanded_coset_domain.ifft(&validity_poly_evals);
    // println!("validity_poly_evals:{:?}", validity_poly_evals);
    // println!("validity_poly_coeffs:{:?} ", validity_poly_coeffs);

    // 8. The DEEP Technique.

    let deep_test_point = Fr::from(93);
    let deep_test_point_forward = deep_test_point.mul(&domain.group_gen_inv());

    let eval = DensePolynomial::from_coefficients_slice(&trace_polys[0]).evaluate(&deep_test_point);
    let deep_poly1_evals = deep_quotient(
        &trace_zk_commitments[0],
        &deep_test_point,
        &eval,
        &expanded_coset_domain,
    );

    let eval0 =
        DensePolynomial::from_coefficients_slice(&trace_polys[1]).evaluate(&deep_test_point);
    let eval1 = DensePolynomial::from_coefficients_slice(&trace_polys[1])
        .evaluate(&deep_test_point_forward);
    let deep_poly2_evals = deep_quotient_two_points(
        &trace_zk_commitments[1],
        &[deep_test_point, deep_test_point_forward],
        &[eval0, eval1],
        &expanded_coset_domain,
    );

    let eval0 =
        DensePolynomial::from_coefficients_slice(&trace_polys[2]).evaluate(&deep_test_point);
    let eval1 = DensePolynomial::from_coefficients_slice(&trace_polys[2])
        .evaluate(&deep_test_point_forward);
    let deep_poly3_evals = deep_quotient_two_points(
        &trace_zk_commitments[2],
        &[deep_test_point, deep_test_point_forward],
        &[eval0, eval1],
        &expanded_coset_domain,
    );

    let eval = DensePolynomial::from_coefficients_slice(&trace_polys[3]).evaluate(&deep_test_point);
    let deep_poly4_evals = deep_quotient(
        &trace_zk_commitments[3],
        &deep_test_point,
        &eval,
        &expanded_coset_domain,
    );

    let eval = DensePolynomial::from_coefficients_slice(&trace_polys[4]).evaluate(&deep_test_point);
    let deep_poly5_evals = deep_quotient(
        &trace_zk_commitments[4],
        &deep_test_point,
        &eval,
        &expanded_coset_domain,
    );

    let eval = DensePolynomial::from_coefficients_slice(&trace_polys[5]).evaluate(&deep_test_point);
    let deep_poly6_evals = deep_quotient(
        &trace_zk_commitments[5],
        &deep_test_point,
        &eval,
        &expanded_coset_domain,
    );

    let eval =
        DensePolynomial::from_coefficients_slice(&validity_poly_coeffs).evaluate(&deep_test_point);
    let deep_poly7_evals = deep_quotient(
        &validity_poly_evals,
        &deep_test_point,
        &eval,
        &expanded_coset_domain,
    );

    // 9. Mixing for FRI.

    let fri_input_poly_evals = mix(
        &[
            deep_poly1_evals,
            deep_poly2_evals,
            deep_poly3_evals,
            deep_poly4_evals,
            deep_poly5_evals,
            deep_poly6_evals,
            deep_poly7_evals,
        ],
        &Fr::from(21),
    );

    println!(
        "fri_input_poly_evals: {:?}",
        expanded_domain.ifft(&fri_input_poly_evals)
    );
}
