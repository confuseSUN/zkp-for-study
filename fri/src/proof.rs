use std::fmt::Display;

use ark_ff::PrimeField;
use merkle_tree::Proof;

#[derive(Default, Debug)]
pub struct FriProof<T: Display> {
    pub merkle_root: Vec<String>,
    pub last_codeword: Vec<T>,
    pub colinearity_tests: Vec<(T, T, T)>,
    pub merkle_auth_paths: Vec<(Proof<T>, Proof<T>, Proof<T>)>,
}

impl<T: PrimeField> FriProof<T> {
    pub fn push_root(&mut self, root: String) {
        self.merkle_root.push(root)
    }

    pub fn push_colinearity_test(&mut self, colinearity_test: (T, T, T)) {
        self.colinearity_tests.push(colinearity_test)
    }

    pub fn push_merkle_auth_paths(&mut self, merkle_auth_path: (Proof<T>, Proof<T>, Proof<T>)) {
        self.merkle_auth_paths.push(merkle_auth_path)
    }
}
