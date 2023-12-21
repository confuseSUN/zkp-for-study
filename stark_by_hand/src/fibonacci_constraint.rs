use ark_ff::PrimeField;
use std::{ops::*, fmt::Formatter};
use std::fmt::{Debug, self};

use crate::column::{Column, Label};

pub struct FibonacciConstraint<F: PrimeField> {
    pub relation_constraint: Column<F>,
    pub first_input_constraint: Column<F>,
    pub second_input_constraint: Column<F>,
    pub termination_constraint: Column<F>,
    pub first_copy_constraint: Column<F>,
    pub second_copy_constraint: Column<F>,
}

impl<F: PrimeField> FibonacciConstraint<F> {
    pub fn construct(trace: &[Column<F>], public_inputs: &[F], label: Label) -> Self {
        assert_eq!(trace.len(), 6);
        assert_eq!(public_inputs.len(), 3);

        let state0 = trace[2].clone().sub(&trace[0]).sub(&trace[1]);
        let state1 = trace[3].clone().add(&trace[4]).add(&trace[5]);
        let relation_constraint = state0.mul(&state1);

        let state0 = trace[0].subtract_scalar(&public_inputs[0]);
        let first_input_constraint = state0.mul(&trace[3]);

        let state0 = trace[1].subtract_scalar(&public_inputs[1]);
        let second_input_constraint = state0.mul(&trace[3]);

        let state0 = trace[2].subtract_scalar(&public_inputs[2]);
        let termination_constraint = state0.mul(&trace[5]);

        let mut first_copy_constraint = Vec::new();
        for i in 0..trace[0].len() {
            let x = match label {
                Label::TraCeColumn => trace[0].get_raw()[i]
                    .sub(trace[1].get_raw()[(i + trace[0].len() - 1) % trace[0].len()]),
                Label::NotTraCeColumn => trace[0].get_raw()[i]
                    .sub(trace[1].get_raw()[(i + trace[0].len() - 4) % trace[0].len()]),
            };
            first_copy_constraint.push(x.mul(trace[4].get_raw()[i]))
        }

        let mut second_copy_constraint = Vec::new();
        for i in 0..trace[0].len() {
            let x = match label {
                Label::TraCeColumn => trace[1].get_raw()[i]
                    .sub(trace[2].get_raw()[(i + trace[0].len() - 1) % trace[0].len()]),
                Label::NotTraCeColumn => trace[1].get_raw()[i]
                    .sub(trace[2].get_raw()[(i + trace[0].len() - 4) % trace[0].len()]),
            };
            second_copy_constraint.push(x.mul(trace[4].get_raw()[i]))
        }

        Self {
            relation_constraint,
            first_input_constraint,
            second_input_constraint,
            termination_constraint,
            first_copy_constraint: Column::from(&first_copy_constraint),
            second_copy_constraint: Column::from(&second_copy_constraint),
        }
    }

    pub fn mix(&self, alpha: &F) -> Column<F> {
        let mut res = self.relation_constraint.clone();

        res.add_assign(&self.first_input_constraint.scale(alpha));

        let alp_pow_2 = alpha.mul(alpha);
        res.add_assign(&&self.second_input_constraint.scale(&alp_pow_2));

        let alp_pow_3 = alp_pow_2.mul(alpha);
        res.add_assign(&&&self.termination_constraint.scale(&alp_pow_3));

        let alp_pow_4 = alp_pow_3.mul(alpha);
        res.add_assign(&&&&self.first_copy_constraint.scale(&alp_pow_4));

        let alp_pow_5 = alp_pow_4.mul(alpha);
        res.add_assign(&&&&&self.second_copy_constraint.scale(&alp_pow_5));

        res
    }
}

impl<F: PrimeField> Debug for FibonacciConstraint<F>{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("relation: [")?;
        for (i, x) in self.relation_constraint.get_raw_ref().iter().enumerate() {
            let v = x.into_bigint().to_string();
            if i != self.relation_constraint.len() - 1 {
                f.write_fmt(format_args!("{},", v))?;
            } else {
                f.write_fmt(format_args!("{}", v))?;
            }
        }
        f.write_str("] \n")?;

        f.write_str("first input: [")?;
        for (i, x) in self.first_input_constraint.get_raw_ref().iter().enumerate() {
            let v = x.into_bigint().to_string();
            if i != self.first_input_constraint.len() - 1 {
                f.write_fmt(format_args!("{},", v))?;
            } else {
                f.write_fmt(format_args!("{}", v))?;
            }
        }
        f.write_str("] \n")?;

        f.write_str("second input: [")?;
        for (i, x) in self.second_input_constraint.get_raw_ref().iter().enumerate() {
            let v = x.into_bigint().to_string();
            if i != self.second_input_constraint.len() - 1 {
                f.write_fmt(format_args!("{},", v))?;
            } else {
                f.write_fmt(format_args!("{}", v))?;
            }
        }
        f.write_str("] \n")?;

        f.write_str("termination: [")?;
        for (i, x) in self.termination_constraint.get_raw_ref().iter().enumerate() {
            let v = x.into_bigint().to_string();
            if i != self.termination_constraint.len() - 1 {
                f.write_fmt(format_args!("{},", v))?;
            } else {
                f.write_fmt(format_args!("{}", v))?;
            }
        }
        f.write_str("] \n")?;

        f.write_str("first copy: [")?;
        for (i, x) in self.first_copy_constraint.get_raw_ref().iter().enumerate() {
            let v = x.into_bigint().to_string();
            if i != self.first_copy_constraint.len() - 1 {
                f.write_fmt(format_args!("{},", v))?;
            } else {
                f.write_fmt(format_args!("{}", v))?;
            }
        }
        f.write_str("] \n")?;

        f.write_str("second copy: [")?;
        for (i, x) in self.second_copy_constraint.get_raw_ref().iter().enumerate() {
            let v = x.into_bigint().to_string();
            if i != self.second_copy_constraint.len() - 1 {
                f.write_fmt(format_args!("{},", v))?;
            } else {
                f.write_fmt(format_args!("{}", v))?;
            }
        }
        f.write_str("] \n")
    }
}