use ark_ff::{batch_inversion, PrimeField};
use ark_poly::{EvaluationDomain, Radix2EvaluationDomain};

pub fn vec_add<F: PrimeField>(a: &[F], b: &[F]) -> Vec<F> {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(x, y)| x.add(y)).collect()
}

pub fn vec_sub<F: PrimeField>(a: &[F], b: &[F]) -> Vec<F> {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(x, y)| x.sub(y)).collect()
}

pub fn vec_mul<F: PrimeField>(a: &[F], b: &[F]) -> Vec<F> {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(x, y)| x.mul(y)).collect()
}

pub fn vec_div<F: PrimeField>(a: &[F], b: &[F]) -> Vec<F> {
    assert_eq!(a.len(), b.len());
    let mut b_invs = b.to_vec();
    batch_inversion(&mut b_invs);
    a.iter().zip(b_invs.iter()).map(|(x, y)| x.mul(y)).collect()
}

pub fn mix<F: PrimeField>(data: &[Vec<F>], alpha: &F) -> Vec<F> {
    let mut res = data[0].clone();
    let mut multiplier = alpha.clone();
    for d in data.iter().skip(1) {
        let scale = d.iter().map(|x| x.mul(&multiplier)).collect::<Vec<_>>();
        res = vec_add(&res, &scale);
        multiplier.mul_assign(alpha);
    }
    res
}

pub fn deep_quotient<F: PrimeField>(
    poly: &[F],
    point: &F,
    eval: &F,
    domain: &Radix2EvaluationDomain<F>,
) -> Vec<F> {
    let numerator = poly.iter().map(|x| x.sub(eval)).collect::<Vec<_>>();
    let denominator: Vec<F> = domain.elements().into_iter().map(|x| x - point).collect();
    vec_div(&numerator, &denominator)

    // let mut poly = self.clone();
    // poly.0[0].sub_assign(eval);
    // let numerator = DensePolynomial::from_coefficients_slice(poly.get_raw_ref());

    // let denominator = DensePolynomial::from_coefficients_slice(&[point.neg(),F::ONE]);
    // let res =  numerator.div(&denominator);

    // Self(res.coeffs)
}
