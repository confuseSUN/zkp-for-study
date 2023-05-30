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
    use crate::linear_combination::{LinearCombination, Sign};
    use ark_std::ops::*;
    use sample_field::BN254Fr;

    use super::Variable;

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
}
