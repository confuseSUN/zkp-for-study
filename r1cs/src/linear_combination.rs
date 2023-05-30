use ark_std::ops::*;

use ark_ff::Field;

/// Variable index
pub type VarIndex = usize;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Sign {
    Positive,
    Negative,
}

impl Neg for Sign {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        match self {
            Sign::Positive => Sign::Negative,
            Sign::Negative => Sign::Positive,
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Variable<F: Field> {
    Constant(F),
    Instance(VarIndex),
    Witness(VarIndex),
}

impl<F: Field, L: Into<LinearCombination<F>>> Add<L> for Variable<F> {
    type Output = LinearCombination<F>;

    #[inline]
    fn add(self, rhs: L) -> Self::Output {
        let l = LinearCombination::from(self);
        let r: LinearCombination<F> = rhs.into();
        l + r
    }
}

impl<F: Field, L: Into<LinearCombination<F>>> Sub<L> for Variable<F> {
    type Output = LinearCombination<F>;

    fn sub(self, rhs: L) -> Self::Output {
        let l = LinearCombination::from(self);
        let r: LinearCombination<F> = rhs.into();
        l - r
    }
}

impl<F: Field> Neg for Variable<F> {
    type Output = LinearCombination<F>;

    #[inline]
    fn neg(self) -> Self::Output {
        -LinearCombination::from(self)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct LinearCombination<F: Field> {
    pub terms: Vec<(Variable<F>, Sign)>,
}

impl<F: Field> LinearCombination<F> {
    pub fn eval(&self, witness: &[F]) -> F {
        let mut res = F::zero();

        for (var, sign) in self.terms.iter() {
            match (var, sign) {
                (Variable::Constant(v), Sign::Positive) => res += v,
                (Variable::Constant(v), Sign::Negative) => res -= v,
                (Variable::Instance(i) | Variable::Witness(i), Sign::Positive) => {
                    res += witness[*i]
                }
                (Variable::Instance(i) | Variable::Witness(i), Sign::Negative) => {
                    res -= witness[*i]
                }
            };
        }

        res
    }
}

impl<F: Field> From<Variable<F>> for LinearCombination<F> {
    fn from(v: Variable<F>) -> LinearCombination<F> {
        LinearCombination {
            terms: vec![(v, Sign::Positive)],
        }
    }
}

impl<'a, F: Field> From<&'a Variable<F>> for LinearCombination<F> {
    fn from(v: &'a Variable<F>) -> LinearCombination<F> {
        LinearCombination {
            terms: vec![(v.clone(), Sign::Positive)],
        }
    }
}

impl<F: Field, L: Into<LinearCombination<F>>> Add<L> for LinearCombination<F> {
    type Output = Self;

    #[inline]
    fn add(mut self, rhs: L) -> Self::Output {
        self.terms.extend_from_slice(&rhs.into().terms);
        self
    }
}

impl<F: Field, L: Into<LinearCombination<F>>> Sub<L> for LinearCombination<F> {
    type Output = Self;

    fn sub(mut self, rhs: L) -> Self::Output {
        let tmp = rhs.into().neg();
        self.terms.extend_from_slice(&tmp.terms);
        self
    }
}

impl<F: Field> Neg for LinearCombination<F> {
    type Output = Self;

    #[inline]
    fn neg(mut self) -> Self::Output {
        for (_, ops) in self.terms.iter_mut() {
            *ops = -*ops
        }

        self
    }
}

#[cfg(test)]
mod tests {
    use super::Variable;
    use crate::linear_combination::Variable::*;
    use crate::linear_combination::{LinearCombination, Sign};
    use ark_ff::One;
    use ark_std::ops::*;
    use sample_field::BN254Fr;

    #[test]
    fn test_lc() {
        let var1 = Variable::<BN254Fr>::Witness(1);
        let var2 = Variable::<BN254Fr>::Witness(2);
        let var3: LinearCombination<BN254Fr> = var1 + var2;
        assert_eq!(
            var3.terms,
            vec![(var1, Sign::Positive), (var2, Sign::Positive)]
        );

        let var1_neg = var1.neg();
        assert_eq!(var1_neg.terms, vec![(var1, Sign::Negative)]);

        let var4: LinearCombination<BN254Fr> = var3.clone() - var1;
        let var5: LinearCombination<BN254Fr> = var3.clone() + var1_neg;

        assert_eq!(var4, var5)
    }

    #[test]
    fn test_eval() {
        // [1, 3, 7, 21, 9, 27, 14, 9]
        let witness = [
            BN254Fr::one(),
            BN254Fr::from(3),
            BN254Fr::from(7),
            BN254Fr::from(21),
            BN254Fr::from(9),
            BN254Fr::from(27),
            BN254Fr::from(14),
            BN254Fr::from(9),
        ];

        let lc = LinearCombination {
            terms: vec![
                (Witness::<BN254Fr>(5), Sign::Positive),
                (Witness::<BN254Fr>(3), Sign::Negative),
                (Witness::<BN254Fr>(1), Sign::Negative),
                (Witness::<BN254Fr>(6), Sign::Positive),
                (Constant(BN254Fr::from(8)), Sign::Negative),
            ],
        };

        let eval = lc.eval(&witness);

        // = 27 - 21 - 3 + 14 -8
        assert_eq!(eval, BN254Fr::from(9))
    }
}
