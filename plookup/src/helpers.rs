use std::ops::{Add, Mul, MulAssign, Sub};

use ark_ff::{Field, One, Zero};
use ark_poly::{EvaluationDomain, Polynomial, Radix2EvaluationDomain, UVPolynomial};
use kzg::{Poly, Scalar};

use crate::table::SampleTable;

pub fn compute_n_lagrange_poly<E: EvaluationDomain<Scalar>>(domain: &E, n: usize) -> Poly {
    let mut evals = vec![Scalar::zero(); domain.size()];
    evals[n] = Scalar::one();
    domain.ifft_in_place(&mut evals);
    Poly::from_coefficients_vec(evals)
}

pub fn compute_t_poly(
    z_poly: &Poly,
    f_poly: &Poly,
    t_poly: &Poly,
    h1_poly: &Poly,
    h2_poly: &Poly,
    domain: &Radix2EvaluationDomain<Scalar>,
    beta: &Scalar,
    gamma: &Scalar,
) -> Poly {
    let domain_4n = Radix2EvaluationDomain::<Scalar>::new(domain.size() * 4).unwrap();
    let group_4n = domain_4n.elements();
    let gn = domain.elements().last().unwrap();
    let l1_poly = compute_n_lagrange_poly(domain, 0);
    let ln_poly = compute_n_lagrange_poly(domain, domain.size() - 1);
    let mut h1_evals = domain_4n.fft(&h1_poly);
    let mut h2_evals = domain_4n.fft(&h2_poly);
    let ln_evals = domain_4n.fft(&ln_poly);
    let mut z_evals = domain_4n.fft(&z_poly);
    let f_evals = domain_4n.fft(&f_poly);
    let mut t_evals = domain_4n.fft(&t_poly);
    for i in 0..4 {
        h1_evals.push(h1_evals[i]);
        h2_evals.push(h2_evals[i]);
        t_evals.push(t_evals[i]);
        z_evals.push(z_evals[i]);
    }

    let term1_poly = l1_poly
        .add(ln_poly.clone())
        .mul(&z_poly.sub(&Poly::from_coefficients_vec(vec![Scalar::one()])));

    let (_, remainder) = term1_poly.divide_by_vanishing_poly(*domain).unwrap();
    assert!(remainder.is_empty());

    let mut term2_evals = vec![];
    for i in 0..domain_4n.size() {
        let eval = ln_evals[i].mul(h1_evals[i].sub(h2_evals[i + 4]));
        term2_evals.push(eval);
    }
    let term2_coefs = domain_4n.ifft(&term2_evals);
    let term2_poly = Poly::from_coefficients_vec(term2_coefs);

    let (_, remainder) = term2_poly.divide_by_vanishing_poly(*domain).unwrap();
    assert!(remainder.is_empty());

    let beta_plus_one = beta.add(Scalar::one());
    let gamma_mul_beta_plus_one = gamma.mul(beta_plus_one);
    let term3_evals = (0..domain_4n.size())
        .into_iter()
        .zip(group_4n)
        .map(|(i, ele)| {
            let tmp = ele.sub(gn);

            let mut term1 = tmp.mul(z_evals[i].mul(beta_plus_one));
            term1.mul_assign(gamma.add(f_evals[i]));
            term1.mul_assign(
                gamma_mul_beta_plus_one
                    .add(t_evals[i])
                    .add(t_evals[i + 4].mul(beta)),
            );

            let mut term2 = tmp.mul(z_evals[i + 4]);
            term2.mul_assign(
                gamma_mul_beta_plus_one
                    .add(h1_evals[i])
                    .add(h1_evals[i + 4].mul(beta)),
            );
            term2.mul_assign(
                gamma_mul_beta_plus_one
                    .add(h2_evals[i])
                    .add(h2_evals[i + 4].mul(beta)),
            );

            term1.sub(term2)
        })
        .collect::<Vec<Scalar>>();

    let term3_coefs = domain_4n.ifft(&term3_evals);
    let term3_poly = Poly::from_coefficients_vec(term3_coefs);

    let (_, remainder) = term3_poly.divide_by_vanishing_poly(*domain).unwrap();
    assert!(remainder.is_empty());

    let term_sum = term1_poly.add(term2_poly).add(term3_poly);
    let (t_poly, remainder) = term_sum.divide_by_vanishing_poly(*domain).unwrap();
    assert!(remainder.is_empty());

    t_poly
}

pub fn compute_r_poly(
    z_poly: &Poly,
    h1_poly: &Poly,
    h2_poly: &Poly,
    quotient_poly: &Poly,
    f_poly_eval_zeta: &Scalar,
    t_poly_eval_zeta: &Scalar,
    t_poly_eval_zeta_omega: &Scalar,
    zeta: &Scalar,
    beta: &Scalar,
    gamma: &Scalar,
    zeta_omega: &Scalar,
    domain: &Radix2EvaluationDomain<Scalar>,
) -> Poly {
    let h1_poly_eval_zeta = h1_poly.evaluate(&zeta);
    let h1_poly_eval_zeta_omega = h1_poly.evaluate(&zeta_omega);
    let h2_poly_eval_zeta_omega = h2_poly.evaluate(&zeta_omega);
    let z_poly_eval_zeta_omega = z_poly.evaluate(&zeta_omega);
    let gn = domain.elements().last().unwrap();
    let beta_plus_one = beta.add(Scalar::one());
    let gamma_mul_beta_plus_one = gamma.mul(beta_plus_one);
    let l1_poly = compute_n_lagrange_poly(domain, 0);
    let ln_poly = compute_n_lagrange_poly(domain, domain.size() - 1);
    let l1_poly_eval_zeta = l1_poly.evaluate(zeta);
    let ln_poly_eval_zeta = ln_poly.evaluate(zeta);

    let term1 = z_poly
        .sub(&Poly::from_coefficients_vec(vec![Scalar::one()]))
        .mul(l1_poly_eval_zeta.add(ln_poly_eval_zeta));

    let term2 = z_poly
        .mul(zeta.sub(gn))
        .mul(beta_plus_one)
        .mul(gamma.add(f_poly_eval_zeta))
        .mul(
            gamma_mul_beta_plus_one
                .add(t_poly_eval_zeta)
                .add(beta.mul(t_poly_eval_zeta_omega)),
        );

    let term3 = h2_poly.add(&Poly::from_coefficients_vec(vec![
        gamma_mul_beta_plus_one.add(beta.mul(h2_poly_eval_zeta_omega))
    ]));
    let term3 = term3.mul(zeta.sub(gn).mul(z_poly_eval_zeta_omega)).mul(
        gamma_mul_beta_plus_one
            .add(h1_poly_eval_zeta)
            .add(beta.mul(h1_poly_eval_zeta_omega)),
    );

    let ln_poly_eval_zeta = ln_poly.evaluate(zeta);
    let term4 = h1_poly
        .sub(&Poly::from_coefficients_vec(vec![h2_poly_eval_zeta_omega]))
        .mul(ln_poly_eval_zeta);

    let r_poly = term1
        .add(term2)
        .sub(&term3)
        .add(term4)
        .sub(&quotient_poly.mul(domain.vanishing_polynomial().evaluate(&zeta)));
    r_poly
}

pub fn compute_z_poly(
    n: usize,
    f_table: &SampleTable,
    t_table: &SampleTable,
    h1_table: &SampleTable,
    h2_table: &SampleTable,
    gamma: &Scalar,
    beta: &Scalar,
    domain: &Radix2EvaluationDomain<Scalar>,
) -> Poly {
    let mut z_evals = vec![Scalar::one()];
    let mut numerator = Scalar::one();
    let mut denominator = Scalar::one();

    let beta_plus_one = beta.add(Scalar::one());
    let gamma_mul_beta_plus_one = gamma.mul(&beta_plus_one);

    for i in 0..n {
        // numerator = (\beta + 1) * (\gamma + f_i) * (\gamma * (\beta + 1) + t_i + \beta * t_{i+1})
        numerator.mul_assign(beta_plus_one);
        numerator.mul_assign(gamma.add(f_table.0[i]));
        numerator.mul_assign(
            gamma_mul_beta_plus_one
                .add(t_table.0[i])
                .add(beta.mul(t_table.0[i + 1])),
        );

        // denominator = (\gamma * (\beta + 1) + h1_i + \beta * h1_{i+1}) * (\gamma * (\beta + 1) + h2_i + \beta * h2_{i+1})
        let term1 = gamma_mul_beta_plus_one
            .add(h1_table.0[i])
            .add(beta.mul(h1_table.0[i + 1]));
        let term2 = gamma_mul_beta_plus_one
            .add(h2_table.0[i])
            .add(beta.mul(h2_table.0[i + 1]));
        denominator.mul_assign(term1.mul(&term2));

        let denominator_inv = denominator.inverse().unwrap();
        z_evals.push(numerator.mul(&denominator_inv));
    }
    assert_eq!(z_evals.last().unwrap().to_owned(), Scalar::one());

    let z_coefs = domain.ifft(&z_evals);
    let z_poly = Poly::from_coefficients_vec(z_coefs);
    z_poly
}
