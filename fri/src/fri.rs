use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial, Polynomial};
use ark_std::rand::Rng;
use merkle_tree::MerkleTree;
use merlin::Transcript;

use crate::{proof::FriProof, transcripts::GlobalTranscript, utils::test_rng_helper_from_scalar};

pub struct FRI<F: PrimeField> {
    pub offset: F,
    pub omega: F,
    pub expansion_factor: usize,
    pub num_colinearity_tests: usize,
    pub codeword_length: usize,
}

impl<F: PrimeField> FRI<F> {
    pub fn new(
        codeword_length: usize,
        expansion_factor: usize,
        num_colinearity_tests: usize,
    ) -> Self {
        let offset = F::GENERATOR;
        let omega = F::get_root_of_unity(codeword_length as u64).unwrap();
        assert!(omega.pow(&[codeword_length as u64]) == F::one());

        FRI {
            offset,
            omega,
            expansion_factor,
            num_colinearity_tests,
            codeword_length,
        }
    }

    pub fn num_rounds(&self) -> usize {
        let mut num_rounds = 0;
        let mut codeword_length = self.codeword_length;

        while codeword_length > self.expansion_factor
            && 4 * self.num_colinearity_tests < codeword_length
        {
            codeword_length /= 2;
            num_rounds += 1;
        }

        num_rounds
    }

    pub fn prove(&self, codeword: &[F]) -> FriProof<F> {
        assert!(codeword.len() == self.codeword_length);

        let mut transcript = self.init_transcript();
        let (mut proof, codewords) = self.commit(codeword, &mut transcript);
        
        let mut indexs = self.sample_index(
            codewords[1].len(),
            codewords[codewords.len() - 1].len(),
            &mut transcript,
        );

        for i in 0..(codewords.len() - 1) {
            indexs = indexs
                .into_iter()
                .map(|x| x % (codewords[i].len() / 2))
                .collect();
            self.query(&codewords[i], &codewords[i + 1], &indexs, &mut proof);
        }

        proof
    }

    pub fn commit(
        &self,
        codeword: &[F],
        transcript: &mut Transcript,
    ) -> (FriProof<F>, Vec<Vec<F>>) {
        let mut omega = self.omega;
        let mut offset = self.offset;
        let mut proof = FriProof::default();
        let mut codeword = codeword.to_vec();
        let mut codewords = vec![];
        let one = F::one();
        let two = one.add(&one);
        let two_inverse = two.inverse().unwrap();

        for r in 0..self.num_rounds() {
            let codeword_length = codeword.len();
            assert!(omega.pow(&[(codeword_length - 1) as u64]) == omega.inverse().unwrap());

            let db = MerkleTree::<F>::from_vec(codeword.to_vec());
            let root = db.root_hash().unwrap();
            proof.push_root(root.to_string());

            transcript.append_message(b"root", root.as_bytes());

            codewords.push(codeword.clone());

            if r != self.num_rounds() - 1 {
                let alpha: F = transcript.get_challenge(b"alpha");

                let mut next_codeword = vec![];

                for i in 0..codeword_length / 2 {
                    let omega_pow_i = omega.pow(&[i as u64]);
                    let omega_pow_i_inv = omega_pow_i.inverse().unwrap();
                    let offset_inv = offset.inverse().unwrap();
                    let temp1 = one
                        .add(alpha.mul(offset_inv.mul(&omega_pow_i_inv)))
                        .mul(&codeword[i]);
                    let temp2 = one
                        .sub(alpha.mul(offset_inv.mul(&omega_pow_i_inv)))
                        .mul(&codeword[i + codeword_length / 2]);

                    next_codeword.push(two_inverse.mul(temp1.add(&temp2)));
                }

                omega = omega.square();
                offset = offset.square();
                codeword = next_codeword;
            }
        }

        transcript.append_scalars(&codeword);
        proof.last_codeword = codeword;

        (proof, codewords)
    }

    pub fn sample_index(
        &self,
        size: usize,
        reduced_size: usize,
        transcript: &mut Transcript,
    ) -> Vec<usize> {
        assert!(2 * reduced_size >= self.num_colinearity_tests);
        let mut num_colinearity_tests = self.num_colinearity_tests;
        let mut indexs = vec![];
        let mut re_indexs = vec![];

        let seed: F = transcript.get_challenge(b"seed");
        let mut rng = test_rng_helper_from_scalar(&seed);

        while num_colinearity_tests > 0 {
            let index = rng.gen_range(0..size);
            let re_index = index % reduced_size;
            if !re_indexs.contains(&re_index) {
                indexs.push(index);
                re_indexs.push(re_index);
                num_colinearity_tests -= 1;
            }
        }

        indexs
    }

    pub fn query(
        &self,
        current_codeword: &[F],
        next_codeword: &[F],
        indexs: &[usize],
        proof: &mut FriProof<F>,
    ) {
        let first_indexs = indexs;
        let second_indexs = indexs
            .into_iter()
            .map(|x| x + (current_codeword.len() / 2))
            .collect::<Vec<usize>>();

        let mut colinearity_tests = vec![];
        let mut merkle_auth_paths = vec![];

        for i in 0..self.num_colinearity_tests {
            // 1. colinearity tests
            let colinearity_test = (
                current_codeword[first_indexs[i]],
                current_codeword[second_indexs[i]],
                next_codeword[first_indexs[i]],
            );
            colinearity_tests.push(colinearity_test);

            // 2. merkle authentication paths
            let current_codeword_db = MerkleTree::<F>::from_vec(current_codeword.to_vec());
            let current_codeword_path1 =
                current_codeword_db.get_proof(current_codeword[first_indexs[i]]);
            let current_codeword_path2 =
                current_codeword_db.get_proof(current_codeword[second_indexs[i]]);

            let next_codeword_db = MerkleTree::<F>::from_vec(next_codeword.to_vec());
            let next_codeword_path = next_codeword_db.get_proof(next_codeword[first_indexs[i]]);
            merkle_auth_paths.push((
                current_codeword_path1,
                current_codeword_path2,
                next_codeword_path,
            ));
        }

        proof.push_colinearity_test(colinearity_tests);
        proof.push_merkle_auth_paths(merkle_auth_paths);
    }

    pub fn verify(&self, proof: &FriProof<F>) {
        let last_codeword_db = MerkleTree::<F>::from_vec(proof.last_codeword.clone());
        assert!(last_codeword_db.root_hash().unwrap() == proof.merkle_root.last().unwrap());

        let mut omega = self.omega;
        let mut offset = self.offset;
        let mut last_omega = omega;
        let mut last_offset = offset;
        for _ in 0..self.num_rounds() - 1 {
            last_omega = last_omega.square();
            last_offset = last_offset.square();
        }
        let last_codeword_length = proof.last_codeword.len();
        assert!(
            last_omega.inverse().unwrap() == last_omega.pow(&[(last_codeword_length - 1) as u64])
        );

        //    let last_domain  = Radix2EvaluationDomain::<T>::new(last_codeword_length).unwrap();
        //    assert!(last_domain.group_gen == last_omega);
        //    let last_poly_coefs = last_domain.coset_ifft(&proof.last_codeword);
        //    let last_poly = DensePolynomial::from_coefficients_vec(last_poly_coefs).mul(last_offset.inverse().unwrap());
        //    let  degree = last_codeword_length / self.expansion_factor - 1;
        //    assert!(last_poly.degree() <= degree )

        let mut alphas = vec![];
        let mut transcript = self.init_transcript();
        for r in 0..self.num_rounds() {
            transcript.append_message(b"root", proof.merkle_root[r].as_bytes());

            if r != self.num_rounds() - 1 {
                let alpha: F = transcript.get_challenge(b"alpha");
                alphas.push(alpha);
            }
        }
        transcript.append_scalars(&proof.last_codeword);

        let top_indexs = self.sample_index(
            self.codeword_length >> 1,
            self.codeword_length >> (self.num_rounds() - 1),
            &mut transcript,
        );

        for r in 0..self.num_rounds() - 1 {
            let first_indexs = top_indexs
                .iter()
                .map(|x| x % (self.codeword_length >> (r + 1)))
                .collect::<Vec<usize>>();
            let second_indexs = first_indexs
                .iter()
                .map(|x| x + (self.codeword_length >> (r + 1)))
                .collect::<Vec<usize>>();

            // 1. verify colinearity test
            for i in 0..self.num_colinearity_tests {
                let x1 = omega.pow(&[first_indexs[i] as u64]).mul(&offset);
                let x2 = omega.pow(&[second_indexs[i] as u64]).mul(&offset);
                let x3 = alphas[r];

                let (y1, y2, y3) = proof.colinearity_tests[r][i];

                Self::colinearity_test((x1, y1), (x2, y2), (x3, y3));
            }

            omega = omega.square();
            offset = offset.square();

            // 2. verify merkle path
            for i in 0..self.num_colinearity_tests {
                let (a, b, c) = &proof.merkle_auth_paths[r][i];
                assert!(a.validate(&proof.merkle_root[r]));
                assert!(b.validate(&proof.merkle_root[r]));
                assert!(c.validate(&proof.merkle_root[r + 1]));
            }
        }
    }

    fn colinearity_test(a: (F, F), b: (F, F), c: (F, F)) {
        let x1_sub_x2 = a.0.sub(&b.0);
        let y1_sub_y2 = a.1.sub(&b.1);
        let k = x1_sub_x2.inverse().unwrap().mul(&y1_sub_y2);
        let q = a.1.sub(&k.mul(a.0));

        let poly = DensePolynomial::from_coefficients_vec(vec![q, k]);
        let eval = poly.evaluate(&c.0);
        assert!(eval == c.1)
    }

    fn init_transcript(&self) -> Transcript {
        let mut transcript = Transcript::new(b"fri");
        transcript.append_u64(b"codeword_length", self.codeword_length as u64);
        transcript.append_u64(b"expansion_factor", self.expansion_factor as u64);
        transcript.append_u64(b"num_colinearity_tests", self.num_colinearity_tests as u64);
        transcript
    }
}
