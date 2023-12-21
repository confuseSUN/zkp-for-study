use ark_ff::FftField;
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain,
};
use column::Label;
use field::Fr;

use crate::column::Column;
use crate::fibonacci_constraint::FibonacciConstraint;

pub mod column;
pub mod fibonacci_constraint;
pub mod field;

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

    println!("The trace polynomial:");
    let mut trace_polys = Vec::new();
    let domain = Radix2EvaluationDomain::<Fr>::new(trace[0].len()).unwrap();
    for values in trace.iter() {
        let poly = Column::from(&domain.ifft(values.get_raw_ref()));
        println!("{:?}", poly);
        trace_polys.push(poly);
    }
    println!();

    println!("The reedsolomon codes of trace polynomial:");
    let mut trace_reedsolomon_codes = Vec::new();
    let expanded_domain = Radix2EvaluationDomain::<Fr>::new(8 * 4).unwrap();
    for coeffs in trace_polys.iter() {
        let poly = DensePolynomial::from_coefficients_slice(coeffs.get_raw_ref());
        let code = Column::from(&poly.evaluate_over_domain(expanded_domain).evals);
        println!("{:?}", code);
        trace_reedsolomon_codes.push(code);
    }
    println!();

    println!("The zk codes of trace polynomial:");
    let mut trace_zk_codes = Vec::with_capacity(trace_polys.len());
    let expanded_coset_domain = expanded_domain.get_coset(Fr::GENERATOR).unwrap();
    for coeffs in trace_polys.iter() {
        let poly = DensePolynomial::from_coefficients_slice(coeffs.get_raw_ref());
        let code = Column::from(&poly.evaluate_over_domain(expanded_coset_domain).evals);
        println!("{:?}", code);
        trace_zk_codes.push(code);
    }
    println!();

    // let x = expanded_coset_domain.ifft(&trace_zk_codes[0]);
    // let p = DensePolynomial::from_coefficients_slice(&x);
    // let r = p.evaluate_over_domain(domain);
    // println("", &r.evals);

    let constraint_column =
        FibonacciConstraint::construct(&trace, &public_inputs, Label::TraCeColumn);
    let constraint_reedsolomon_codes = FibonacciConstraint::construct(
        &trace_reedsolomon_codes,
        &public_inputs,
        Label::NotTraCeColumn,
    );
    let constraint_zk_codes =
        FibonacciConstraint::construct(&trace_zk_codes, &public_inputs, Label::NotTraCeColumn);
    println!("constraint column:");
    println!("equal constraint:{:?}", constraint_column.equal_constraint);
    println!(
        "first input constraint:{:?}",
        constraint_column.first_input_constraint
    );
    println!(
        "second input constraint:{:?}",
        constraint_column.second_input_constraint
    );
    println!(
        "termination input constraint:{:?}",
        constraint_column.termination_constraint
    );
    println!(
        "first permutation constraint:{:?} ",
        constraint_column.first_permutation_constraint
    );
    println!(
        "second permutation constraint:{:?} \n",
        constraint_column.second_permutation_constraint
    );

    println!("constraint reedsolomon codes:");
    println!(
        "equal constraint:{:?}",
        constraint_reedsolomon_codes.equal_constraint
    );
    println!(
        "first input constraint:{:?}",
        constraint_reedsolomon_codes.first_input_constraint
    );
    println!(
        "second input constraint:{:?}",
        constraint_reedsolomon_codes.second_input_constraint
    );
    println!(
        "termination input constraint:{:?}",
        constraint_reedsolomon_codes.termination_constraint
    );
    println!(
        "first permutation constraint:{:?} ",
        constraint_reedsolomon_codes.first_permutation_constraint
    );
    println!(
        "second permutation constraint:{:?} \n",
        constraint_reedsolomon_codes.second_permutation_constraint
    );

    println!("constraint zk codes:");
    println!(
        "equal constraint:{:?}",
        constraint_zk_codes.equal_constraint
    );
    println!(
        "first input constraint:{:?}",
        constraint_zk_codes.first_input_constraint
    );
    println!(
        "second input constraint:{:?}",
        constraint_zk_codes.second_input_constraint
    );
    println!(
        "termination input constraint:{:?}",
        constraint_zk_codes.termination_constraint
    );
    println!(
        "first permutation constraint:{:?} ",
        constraint_zk_codes.first_permutation_constraint
    );
    println!(
        "second permutation constraint:{:?} \n",
        constraint_zk_codes.second_permutation_constraint
    );

    let mix_alpha = Fr::from(3);
    let mixed_constraint_columns = constraint_column.mix(&mix_alpha);
    let mixed_constraint_reedsolomon_codes = constraint_reedsolomon_codes.mix(&mix_alpha);
    let mixed_constraint_zk_codes = constraint_zk_codes.mix(&mix_alpha);

    let equal_constraint_polys = Column::from(
        &expanded_domain.ifft(constraint_reedsolomon_codes.equal_constraint.get_raw_ref()),
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
    let first_permutation_constraint_polys = Column::from(
        &expanded_domain.ifft(
            constraint_reedsolomon_codes
                .first_permutation_constraint
                .get_raw_ref(),
        ),
    );
    let second_permutation_polys = Column::from(
        &expanded_domain.ifft(
            constraint_reedsolomon_codes
                .second_permutation_constraint
                .get_raw_ref(),
        ),
    );
    let mixed_constraint_polys =
        Column::from(&expanded_domain.ifft(mixed_constraint_reedsolomon_codes.get_raw_ref()));
    println!("constraint polys:");
    println!("equal constraint:{:?}", equal_constraint_polys);
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
        first_permutation_constraint_polys
    );
    println!(
        "second permutation constraint:{:?}\n",
        second_permutation_polys
    );

    println!("mixed constraint:{:?} \n", mixed_constraint_polys);

    let z_x = domain.vanishing_polynomial();
    let z_x_reedsolomon_codes =
        Column::from(&z_x.clone().evaluate_over_domain(expanded_domain).evals);
    let z_x_zk_codes = Column::from(&z_x.evaluate_over_domain(expanded_coset_domain).evals);
    println!("zx reedsolomon codes:{:?}", z_x_reedsolomon_codes);
    println!("zx zk codes:{:?} \n", z_x_zk_codes);
}
