use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;

pub trait PolynomialOp<F: PrimeField> {
    fn from_u64(coeffs: &[u64]) -> Self;

    fn subtract_scalar(&self, a: &F) -> Self;
}

impl<F: PrimeField> PolynomialOp<F> for DensePolynomial<F> {
    fn from_u64(coeffs: &[u64]) -> Self {
        let coeffs = coeffs.iter().map(|x| F::from(*x)).collect();
        Self::from_coefficients_vec(coeffs)
    }

    fn subtract_scalar(&self, a: &F) -> Self {
        let mut poly = self.clone();
        poly.iter_mut().for_each(|x| x.sub_assign(a));
        poly
    }
}
