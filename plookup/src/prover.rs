use ark_bls12_381::Fr;
use ark_ff::Zero;
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Polynomial,
    Radix2EvaluationDomain,
};
use kzg::commitment::KZGCommitmentScheme;
use merlin::Transcript;

use crate::{
    helpers::{compute_r_poly, compute_t_poly, compute_z_poly},
    table::{PreProcessedTable, SampleTable},
    transcripts::GlobalTranscript,
    verifier::PlookUpProof,
};
use std::ops::Mul;

pub(super) fn prove(
    mut f_table: SampleTable,
    t_preprocess_table: PreProcessedTable,
    domain: &Radix2EvaluationDomain<Fr>,
    kzg_comm_scheme: &KZGCommitmentScheme,
) -> PlookUpProof {
    let mut transcript = Transcript::new(b"plookup");
    transcript.append_u64(b"size", domain.size as u64);

    let (t_poly, t_comm, t_table) = (
        t_preprocess_table.poly,
        t_preprocess_table.comm,
        t_preprocess_table.table,
    );
    transcript.append_commitent(&t_comm);

    // 1.Compute f_table sorted by t_table.
    let sorted_table = f_table.sort_by(&t_table);

    let f_coefs = domain.ifft(&f_table.0);
    let f_poly = DensePolynomial::from_coefficients_vec(f_coefs);
    let f_comm = kzg_comm_scheme.commit(&f_poly);
    transcript.append_commitent(&f_comm);

    // 2.Compute two helper polynomial(h1, h2) and commit them.
    let n = sorted_table.size() / 2;
    let h1_table = SampleTable::from_scalar(sorted_table.0[..=n].to_vec());
    let h1_coefs = domain.ifft(&h1_table.0);
    let h1_poly = DensePolynomial::from_coefficients_vec(h1_coefs);
    let h1_comm = kzg_comm_scheme.commit(&h1_poly);
    transcript.append_commitent(&h1_comm);

    let h2_table = SampleTable::from_scalar(sorted_table.0[n..].to_vec());
    let h2_coefs = domain.ifft(&h2_table.0);
    let h2_poly = DensePolynomial::from_coefficients_vec(h2_coefs);
    let h2_comm = kzg_comm_scheme.commit(&h2_poly);
    transcript.append_commitent(&h2_comm);

    // 3.Get challenge beta and gamma.
    let beta = transcript.get_challenge(b"beta");
    let gamma = transcript.get_challenge(b"gamma");

    // 4.Compute z(X) and commit it.
    let z_poly = compute_z_poly(
        n, &f_table, &t_table, &h1_table, &h2_table, &gamma, &beta, domain,
    );
    let z_comm = kzg_comm_scheme.commit(&z_poly);
    transcript.append_commitent(&z_comm);

    // 5.Compute quotient polynomial(t(x)).
    let quotient_poly = compute_t_poly(
        &z_poly, &f_poly, &t_poly, &h1_poly, &h2_poly, domain, &beta, &gamma,
    );

    let zeta = transcript.get_challenge(b"zeta");
    let zeta_omega = zeta.mul(&domain.group_gen);
    let f_poly_eval_zeta = f_poly.evaluate(&zeta);
    let t_poly_eval_zeta = t_poly.evaluate(&zeta);
    let t_poly_eval_zeta_omega = t_poly.evaluate(&zeta_omega);

    // 6.Compute r(X) and commit it.
    let r_poly = compute_r_poly(
        &z_poly,
        &h1_poly,
        &h2_poly,
        &quotient_poly,
        &f_poly_eval_zeta,
        &t_poly_eval_zeta,
        &t_poly_eval_zeta_omega,
        &zeta,
        &beta,
        &gamma,
        &zeta_omega,
        domain,
    );
    assert!(r_poly.evaluate(&zeta) == Fr::zero());
    let r_comm = kzg_comm_scheme.commit(&r_poly);
    transcript.append_commitent(&r_comm);

    // 7.Compute two opening proof polynomial.
    let v = transcript.get_challenge(b"v");
    let w_zeta = kzg_comm_scheme.batch_prove(
        &[f_poly, t_poly.clone(), h1_poly.clone(), r_poly],
        &zeta,
        &v,
    );

    let w_zeta_omega =
        kzg_comm_scheme.batch_prove(&[t_poly, h1_poly, h2_poly, z_poly], &zeta_omega, &v);

    PlookUpProof {
        f_comm,
        t_comm,
        h1_comm,
        h2_comm,
        r_comm,
        z_comm,
        opening_witness_zeta: w_zeta,
        opening_witness_zeta_omega: w_zeta_omega,
        domain: domain.to_owned(),
    }
}
