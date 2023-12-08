use ark_ff::PrimeField;

use crate::fft::radix2_fft;

pub fn radix2_ifft<F: PrimeField>(values: &[F], domain: &[F]) -> Vec<F> {
    let mut coeffs = radix2_fft(values, domain);

    let n_inv = F::from(values.len() as u64).inverse().unwrap();
    coeffs.iter_mut().for_each(|x| x.mul_assign(&n_inv));

    coeffs
}
