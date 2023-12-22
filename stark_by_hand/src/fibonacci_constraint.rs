use ark_ff::PrimeField;
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};
use std::fmt::Formatter;
use std::fmt::{self, Debug};

use crate::utils::{mix, vec_add, vec_mul, vec_sub};

pub enum Label {
    Trace,
    ReedsolomonOrZkCommitment,
}

pub struct FibonacciConstraint<F: PrimeField> {
    pub relation_constraint: Vec<F>,
    pub first_input_constraint: Vec<F>,
    pub second_input_constraint: Vec<F>,
    pub termination_constraint: Vec<F>,
    pub first_copy_constraint: Vec<F>,
    pub second_copy_constraint: Vec<F>,
}

impl<F: PrimeField> FibonacciConstraint<F> {
    pub fn construct(trace: &[Vec<F>], public_inputs: &[F], label: Label) -> Self {
        assert_eq!(trace.len(), 6);
        assert_eq!(public_inputs.len(), 3);

        let state0 = vec_sub(&trace[2], &vec_add(&trace[0], &trace[1]));
        let state1 = vec_add(&trace[3], &vec_add(&trace[4], &trace[5]));
        let relation_constraint = vec_mul(&state0, &state1);

        let state = trace[0]
            .iter()
            .map(|x| x.sub(public_inputs[0]))
            .collect::<Vec<F>>();
        let first_input_constraint = vec_mul(&state, &trace[3]);

        let state = trace[1]
            .iter()
            .map(|x| x.sub(public_inputs[1]))
            .collect::<Vec<F>>();
        let second_input_constraint = vec_mul(&state, &trace[3]);

        let state = trace[2]
            .iter()
            .map(|x| x.sub(public_inputs[2]))
            .collect::<Vec<F>>();
        let termination_constraint = vec_mul(&state, &trace[5]);

        let mut first_copy_constraint = Vec::new();
        for i in 0..trace[0].len() {
            let x = match label {
                Label::Trace => {
                    trace[0][i].sub(trace[1][(i + trace[0].len() - 1) % trace[0].len()])
                }
                Label::ReedsolomonOrZkCommitment => {
                    trace[0][i].sub(trace[1][(i + trace[0].len() - 4) % trace[0].len()])
                }
            };
            first_copy_constraint.push(x.mul(trace[4][i]))
        }

        let mut second_copy_constraint = Vec::new();
        for i in 0..trace[0].len() {
            let x = match label {
                Label::Trace => {
                    trace[1][i].sub(trace[2][(i + trace[0].len() - 1) % trace[0].len()])
                }
                Label::ReedsolomonOrZkCommitment => {
                    trace[1][i].sub(trace[2][(i + trace[0].len() - 4) % trace[0].len()])
                }
            };
            second_copy_constraint.push(x.mul(trace[4][i]))
        }

        Self {
            relation_constraint,
            first_input_constraint,
            second_input_constraint,
            termination_constraint,
            first_copy_constraint,
            second_copy_constraint,
        }
    }

    pub fn mix(&self, alpha: &F) -> Vec<F> {
        mix(
            &[
                self.relation_constraint.clone(),
                self.first_input_constraint.clone(),
                self.second_input_constraint.clone(),
                self.termination_constraint.clone(),
                self.first_copy_constraint.clone(),
                self.second_copy_constraint.clone(),
            ],
            alpha,
        )
    }

    pub fn ifft(&self, domain: &Radix2EvaluationDomain<F>) -> Self {
        let relation_constraint = domain.ifft(&self.relation_constraint);
        let first_input_constraint = domain.ifft(&self.first_input_constraint);
        let second_input_constraint = domain.ifft(&self.second_input_constraint);
        let termination_constraint = domain.ifft(&self.termination_constraint);
        let first_copy_constraint = domain.ifft(&self.first_input_constraint);
        let second_copy_constraint = domain.ifft(&self.second_copy_constraint);

        Self {
            relation_constraint,
            first_input_constraint,
            second_input_constraint,
            termination_constraint,
            first_copy_constraint,
            second_copy_constraint,
        }
    }
}

impl<F: PrimeField> Debug for FibonacciConstraint<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("relation: [")?;
        for (i, x) in self.relation_constraint.iter().enumerate() {
            let v = x.into_bigint().to_string();
            if i != self.relation_constraint.len() - 1 {
                f.write_fmt(format_args!("{},", v))?;
            } else {
                f.write_fmt(format_args!("{}", v))?;
            }
        }
        f.write_str("] \n")?;

        f.write_str("first input: [")?;
        for (i, x) in self.first_input_constraint.iter().enumerate() {
            let v = x.into_bigint().to_string();
            if i != self.first_input_constraint.len() - 1 {
                f.write_fmt(format_args!("{},", v))?;
            } else {
                f.write_fmt(format_args!("{}", v))?;
            }
        }
        f.write_str("] \n")?;

        f.write_str("second input: [")?;
        for (i, x) in self.second_input_constraint.iter().enumerate() {
            let v = x.into_bigint().to_string();
            if i != self.second_input_constraint.len() - 1 {
                f.write_fmt(format_args!("{},", v))?;
            } else {
                f.write_fmt(format_args!("{}", v))?;
            }
        }
        f.write_str("] \n")?;

        f.write_str("termination: [")?;
        for (i, x) in self.termination_constraint.iter().enumerate() {
            let v = x.into_bigint().to_string();
            if i != self.termination_constraint.len() - 1 {
                f.write_fmt(format_args!("{},", v))?;
            } else {
                f.write_fmt(format_args!("{}", v))?;
            }
        }
        f.write_str("] \n")?;

        f.write_str("first copy: [")?;
        for (i, x) in self.first_copy_constraint.iter().enumerate() {
            let v = x.into_bigint().to_string();
            if i != self.first_copy_constraint.len() - 1 {
                f.write_fmt(format_args!("{},", v))?;
            } else {
                f.write_fmt(format_args!("{}", v))?;
            }
        }
        f.write_str("] \n")?;

        f.write_str("second copy: [")?;
        for (i, x) in self.second_copy_constraint.iter().enumerate() {
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
