use ark_ff::PrimeField;
use std::ops::*;

use crate::column::{Column, Label};

pub struct FibonacciConstraint<F: PrimeField> {
    pub equal_constraint: Column<F>,
    pub first_input_constraint: Column<F>,
    pub second_input_constraint: Column<F>,
    pub termination_constraint: Column<F>,
    pub first_permutation_constraint: Column<F>,
    pub second_permutation_constraint: Column<F>,
}

impl<F: PrimeField> FibonacciConstraint<F> {
    pub fn construct(trace: &[Column<F>], public_inputs: &[F], label: Label) -> Self {
        assert_eq!(trace.len(), 6);
        assert_eq!(public_inputs.len(), 3);

        let state0 = trace[2].clone().sub(&trace[0]).sub(&trace[1]);
        let state1 = trace[3].clone().add(&trace[4]).add(&trace[5]);
        let equal_constraint = state0.mul(&state1);

        let state0 = trace[0].subtract_scalar(&public_inputs[0]);
        let first_input_constraint = state0.mul(&trace[3]);

        let state0 = trace[1].subtract_scalar(&public_inputs[1]);
        let second_input_constraint = state0.mul(&trace[3]);

        let state0 = trace[2].subtract_scalar(&public_inputs[2]);
        let termination_constraint = state0.mul(&trace[5]);

        let mut first_permutation_constraint = Vec::new();
        for i in 0..trace[0].len() {
            let x = match label {
                Label::TraCeColumn => trace[0].get_raw()[i]
                    .sub(trace[1].get_raw()[(i + trace[0].len() - 1) % trace[0].len()]),
                Label::NotTraCeColumn => trace[0].get_raw()[i]
                    .sub(trace[1].get_raw()[(i + trace[0].len() - 4) % trace[0].len()]),
            };
            first_permutation_constraint.push(x.mul(trace[4].get_raw()[i]))
        }

        let mut second_permutation_constraint = Vec::new();
        for i in 0..trace[0].len() {
            let x = match label {
                Label::TraCeColumn => trace[1].get_raw()[i]
                    .sub(trace[2].get_raw()[(i + trace[0].len() - 1) % trace[0].len()]),
                Label::NotTraCeColumn => trace[1].get_raw()[i]
                    .sub(trace[2].get_raw()[(i + trace[0].len() - 4) % trace[0].len()]),
            };
            second_permutation_constraint.push(x.mul(trace[4].get_raw()[i]))
        }

        Self {
            equal_constraint,
            first_input_constraint,
            second_input_constraint,
            termination_constraint,
            first_permutation_constraint: Column::from(&first_permutation_constraint),
            second_permutation_constraint: Column::from(&second_permutation_constraint),
        }
    }

    pub fn mix(&self, alpha: &F) -> Column<F> {
        let mut res = self.equal_constraint.clone();

        res.add_assign(&self.first_input_constraint.scale(alpha));

        let alp_pow_2 = alpha.mul(alpha);
        res.add_assign(&&self.second_input_constraint.scale(&alp_pow_2));

        let alp_pow_3 = alp_pow_2.mul(alpha);
        res.add_assign(&&&self.termination_constraint.scale(&alp_pow_3));

        let alp_pow_4 = alp_pow_3.mul(alpha);
        res.add_assign(&&&&self.first_permutation_constraint.scale(&alp_pow_4));

        let alp_pow_5 = alp_pow_4.mul(alpha);
        res.add_assign(&&&&&self.second_permutation_constraint.scale(&alp_pow_5));

        res
    }
}
