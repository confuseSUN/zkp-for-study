use ark_ff::PrimeField;
use std::{
    fmt::{self, Debug, Formatter},
    ops::{Add, AddAssign, Mul, Sub},
};

#[derive(Clone)]
pub struct Column<F: PrimeField>(Vec<F>);

impl<F: PrimeField> Column<F> {
    pub fn from(a: &[F]) -> Self {
        Self(a.to_vec())
    }

    pub fn from_u64(a: &[u64]) -> Self {
        let f = a.iter().map(|x| F::from(*x)).collect();
        Self(f)
    }

    pub fn get_raw(&self) -> Vec<F> {
        self.0.clone()
    }

    pub fn get_raw_ref(&self) -> &[F] {
        &self.0
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn subtract_scalar(&self, a: &F) -> Self {
        let mut c = vec![];
        for x in self.0.iter() {
            c.push(x.sub(a))
        }
        Self(c)
    }

    pub fn scale(&self, a: &F) -> Self {
        let mut c = vec![];
        for x in self.0.iter() {
            c.push(x.mul(a))
        }
        Self(c)
    }
}

impl<'a, F: PrimeField> Add<&'a Column<F>> for Column<F> {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self::Output {
        assert_eq!(self.len(), rhs.len());
        let mut r = vec![];
        for (x, y) in self.0.iter().zip(rhs.0.iter()) {
            r.push(x.add(y))
        }
        Self(r)
    }
}

impl<'a, F: PrimeField> AddAssign<&'a Column<F>> for Column<F> {
    fn add_assign(&mut self, rhs: &Self) {
        assert_eq!(self.len(), rhs.len());
        for (x, y) in self.0.iter_mut().zip(rhs.0.iter()) {
            x.add_assign(y)
        }
    }
}

impl<'a, F: PrimeField> Sub<&'a Column<F>> for Column<F> {
    type Output = Self;

    fn sub(self, rhs: &Self) -> Self::Output {
        assert_eq!(self.len(), rhs.len());
        let mut r = vec![];
        for (x, y) in self.0.iter().zip(rhs.0.iter()) {
            r.push(x.sub(y))
        }
        Self(r)
    }
}

impl<'a, F: PrimeField> Mul<&'a Column<F>> for Column<F> {
    type Output = Self;

    fn mul(self, rhs: &Self) -> Self::Output {
        assert_eq!(self.len(), rhs.len());
        let mut r = vec![];
        for (x, y) in self.0.iter().zip(rhs.0.iter()) {
            r.push(x.mul(y))
        }
        Self(r)
    }
}

impl<F: PrimeField> Debug for Column<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("[")?;
        for (i, x) in self.0.iter().enumerate() {
            let v = x.into_bigint().to_string();
            if i != self.len() - 1 {
                f.write_fmt(format_args!("{},", v))?;
            } else {
                f.write_fmt(format_args!("{}", v))?;
            }
        }
        f.write_str("]")
    }
}

pub enum Label {
    TraCeColumn,
    NotTraCeColumn,
}
