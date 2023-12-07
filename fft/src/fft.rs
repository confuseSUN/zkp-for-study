use ark_ff::PrimeField;

use crate::bit_reverse::butterfly;

pub fn radix2_fft<F: PrimeField>(coeffs: &[F], domain: &[F]) -> Vec<F> {
    let n = coeffs.len();
    assert_eq!(n, domain.len() * 2);

    let mut coeffs = coeffs.to_vec();
    butterfly(&mut coeffs);

    // The stage of butterfly operations in FFT.
    let mut stage = 0;

    while (1 << stage) < n {
        let num_domain_elements_needed = 1 << stage;
        let step = num_domain_elements_needed * 2;

        for i in 0..num_domain_elements_needed {
            let w = domain[i * (n / step)];

            for j in (i..n).step_by(step) {
                let w_j = w * coeffs[j + num_domain_elements_needed];

                let r_0 = coeffs[j] + w_j;
                let r_1 = coeffs[j] - w_j;

                coeffs[j] = r_0;
                coeffs[j + num_domain_elements_needed] = r_1;
            }
        }

        stage += 1;
    }

    coeffs
}
