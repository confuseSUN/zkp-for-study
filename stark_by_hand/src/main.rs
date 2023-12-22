use ark_ff::FftField;
use ark_poly::Polynomial;
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain,
};
use column::Label;
use field::Fr;
use std::ops::*;

use crate::column::Column;
use crate::fibonacci_constraint::FibonacciConstraint;
use crate::utils::mix;

pub mod column;
pub mod fibonacci_constraint;
pub mod field;
pub mod utils;
pub mod fri;

fn main() {
    let public_inputs = vec![Fr::from(24), Fr::from(30), Fr::from(28)];

    let data_column_1 = Column::<Fr>::from_u64(&[24, 30, 54, 84, 78, 15, 29, 50]);
    let data_column_2 = Column::<Fr>::from_u64(&[30, 54, 84, 41, 2, 77, 21, 36]);
    let data_column_3 = Column::<Fr>::from_u64(&[54, 84, 41, 28, 71, 17, 92, 33]);

    let control_column_1 = Column::<Fr>::from_u64(&[1, 0, 0, 0, 0, 0, 0, 0]);
    let control_column_2 = Column::<Fr>::from_u64(&[0, 1, 1, 1, 0, 0, 0, 0]);
    let control_column_3 = Column::<Fr>::from_u64(&[0, 0, 0, 1, 0, 0, 0, 0]);

    let trace: Vec<Column<Fr>> = vec![
        data_column_1,
        data_column_2,
        data_column_3,
        control_column_1,
        control_column_2,
        control_column_3,
    ];

    println!("The six trace polynomials:");
    let mut trace_polys = Vec::new();
    let domain = Radix2EvaluationDomain::<Fr>::new(trace[0].len()).unwrap();
    for values in trace.iter() {
        let poly = Column::from(&domain.ifft(values.get_raw_ref()));
        println!("{:?}", poly);
        trace_polys.push(poly);
    }

    println!("\n The reedsolomon codes of trace polynomial:");
    let mut trace_reedsolomon_codes = Vec::new();
    let expanded_domain = Radix2EvaluationDomain::<Fr>::new(8 * 4).unwrap();
    for coeffs in trace_polys.iter() {
        let poly = DensePolynomial::from_coefficients_slice(coeffs.get_raw_ref());
        let evals = Column::from(&poly.evaluate_over_domain(expanded_domain).evals);
        println!("{:?}", evals);
        trace_reedsolomon_codes.push(evals);
    }

    println!("\n The zk commitments of trace polynomial:");
    let mut trace_zk_commitments = Vec::with_capacity(trace_polys.len());
    let expanded_coset_domain = expanded_domain.get_coset(Fr::GENERATOR).unwrap();
    for coeffs in trace_polys.iter() {
        let poly = DensePolynomial::from_coefficients_slice(coeffs.get_raw_ref());
        let evals = Column::from(&poly.evaluate_over_domain(expanded_coset_domain).evals);
        println!("{:?}", evals);
        trace_zk_commitments.push(evals);
    }

    // let x = expanded_coset_domain.ifft(&trace_zk_codes[0]);
    // let p = DensePolynomial::from_coefficients_slice(&x);
    // let r = p.evaluate_over_domain(domain);
    // println("", &r.evals);

    let constraints = FibonacciConstraint::construct(&trace, &public_inputs, Label::TraCeColumn);
    let constraint_reedsolomon_codes = FibonacciConstraint::construct(
        &trace_reedsolomon_codes,
        &public_inputs,
        Label::NotTraCeColumn,
    );
    let constraint_zk_commitments = FibonacciConstraint::construct(
        &trace_zk_commitments,
        &public_inputs,
        Label::NotTraCeColumn,
    );
    println!("\n constraints: \n{:?}", constraints);
    println!("constraint reedsolomon codes:\n{:?}", constraints);
    println!("constraint zk commitments:\n{:?}", constraints);

    let mix_alpha = Fr::from(3);
    let mixed_constraint_column = constraints.mix(&mix_alpha);
    let mixed_constraint_reedsolomon_code = constraint_reedsolomon_codes.mix(&mix_alpha);
    let mixed_constraint_zk_commitment = constraint_zk_commitments.mix(&mix_alpha);

    let relation_constraint_polys = Column::from(
        &expanded_domain.ifft(
            constraint_reedsolomon_codes
                .relation_constraint
                .get_raw_ref(),
        ),
    );
    let first_input_constraint_polys = Column::from(
        &expanded_domain.ifft(
            constraint_reedsolomon_codes
                .first_input_constraint
                .get_raw_ref(),
        ),
    );
    let second_input_constraint_polys = Column::from(
        &expanded_domain.ifft(
            constraint_reedsolomon_codes
                .second_input_constraint
                .get_raw_ref(),
        ),
    );
    let termination_constraint_polys = Column::from(
        &expanded_domain.ifft(
            constraint_reedsolomon_codes
                .termination_constraint
                .get_raw_ref(),
        ),
    );
    let first_copy_constraint_polys = Column::from(
        &expanded_domain.ifft(
            constraint_reedsolomon_codes
                .first_copy_constraint
                .get_raw_ref(),
        ),
    );
    let second_permutation_polys = Column::from(
        &expanded_domain.ifft(
            constraint_reedsolomon_codes
                .second_copy_constraint
                .get_raw_ref(),
        ),
    );
    let mixed_constraint_polys =
        Column::from(&expanded_domain.ifft(mixed_constraint_reedsolomon_code.get_raw_ref()));
    println!("constraint polys:");
    println!("equal constraint:{:?}", relation_constraint_polys);
    println!("first input constraint:{:?}", first_input_constraint_polys);
    println!(
        "second input constraint:{:?}",
        second_input_constraint_polys
    );
    println!(
        "termination input constraint:{:?}",
        termination_constraint_polys
    );
    println!(
        "first permutation constraint:{:?} ",
        first_copy_constraint_polys
    );
    println!(
        "second permutation constraint:{:?}\n",
        second_permutation_polys
    );

    println!("mixed constraint:{:?} \n", mixed_constraint_polys);

    let z_x = domain.vanishing_polynomial();
    let z_x_reedsolomon_code =
        Column::from(&z_x.clone().evaluate_over_domain(expanded_domain).evals);
    let z_x_zk_commitment = Column::from(&z_x.evaluate_over_domain(expanded_coset_domain).evals);
    println!("zx reedsolomon code:{:?}", z_x_reedsolomon_code);
    println!("zx zk commitment:{:?} \n", z_x_zk_commitment);

    // let mut  z_x_zk_commitment_inner = z_x_zk_commitment.get_raw();
    //  batch_inversion(&mut z_x_zk_commitment_inner);
    //  let z_x_zk_commitment_inv = Column::from(&z_x_zk_commitment_inner);
    //  println!("the inversion of zx zk commitment:{:?} \n", z_x_zk_commitment_inv);

    let validity_poly_evals = mixed_constraint_zk_commitment.div(&z_x_zk_commitment);

    println!("validity_poly_evals:{:?} \n", validity_poly_evals);

    let validity_poly_coeffs =
        Column::from(&expanded_coset_domain.ifft(validity_poly_evals.get_raw_ref()));
    println!("validity_poly_coeffs:{:?} \n", validity_poly_coeffs);

    let deep_test_point = Fr::from(93);

    let eval = DensePolynomial::from_coefficients_slice(trace_polys[0].get_raw_ref())
        .evaluate(&deep_test_point);
    let deep_poly1_evals =
        trace_zk_commitments[0].deep_quotient(&deep_test_point, &eval, &expanded_coset_domain);
    let deep_poly2_evals = Column::from_u64(&[
        19, 74, 82, 13, 48, 30, 17, 58, 93, 63, 64, 42, 45, 10, 21, 71, 1, 52, 72, 71, 60, 10, 8,
        44, 41, 5, 14, 16, 49, 15, 78, 41,
    ]);
    let deep_poly3_evals = Column::from_u64(&[
        84, 63, 89, 93, 75, 36, 11, 80, 44, 23, 40, 11, 80, 33, 38, 50, 72, 33, 35, 78, 8, 8, 86,
        36, 9, 25, 88, 95, 65, 22, 50, 91,
    ]);

    let eval = DensePolynomial::from_coefficients_slice(trace_polys[3].get_raw_ref())
        .evaluate(&deep_test_point);
    let deep_poly4_evals =
        trace_zk_commitments[3].deep_quotient(&deep_test_point, &eval, &expanded_coset_domain);

    let eval = DensePolynomial::from_coefficients_slice(trace_polys[4].get_raw_ref())
        .evaluate(&deep_test_point);
    let deep_poly5_evals =
        trace_zk_commitments[4].deep_quotient(&deep_test_point, &eval, &expanded_coset_domain);

    let eval = DensePolynomial::from_coefficients_slice(trace_polys[5].get_raw_ref())
        .evaluate(&deep_test_point);
    let deep_poly6_evals =
        trace_zk_commitments[5].deep_quotient(&deep_test_point, &eval, &expanded_coset_domain);

        let eval = DensePolynomial::from_coefficients_slice(validity_poly_coeffs.get_raw_ref())
        .evaluate(&deep_test_point);
    let deep_poly7_evals =
        validity_poly_evals.deep_quotient(&deep_test_point, &eval, &expanded_coset_domain);

        println!("{:?}",expanded_domain.ifft(deep_poly7_evals.get_raw_ref()));

        let fri_input = mix(&[
            deep_poly1_evals,deep_poly2_evals,deep_poly3_evals,deep_poly4_evals,deep_poly5_evals,deep_poly6_evals,deep_poly7_evals
        ], &Fr::from(21));

        println!("{:?}",expanded_domain.ifft(fri_input.get_raw_ref()));

}
