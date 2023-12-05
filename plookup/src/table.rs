use ark_bls12_381::{Fr, G1Projective};
use ark_poly::{
    univariate::DensePolynomial, DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain,
};
use kzg::commitment::KZGCommitmentScheme;

use crate::{prover, verifier::PlookUpProof};

#[derive(Debug, Clone)]
pub struct SampleTable(pub Vec<Fr>);

impl SampleTable {
    pub fn from_scalar(t: Vec<Fr>) -> Self {
        Self(t)
    }

    pub fn from_u64(t: Vec<u64>) -> Self {
        let t = t.into_iter().map(|x| Fr::from(x)).collect();
        Self(t)
    }

    pub fn pad(&mut self) {
        let size = self.size();
        let power_of_two = size.next_power_of_two();
        let diff = power_of_two - size;

        let ele = vec![*self.0.last().unwrap(); diff];
        self.0.extend(ele);

        assert!(self.size().is_power_of_two())
    }

    pub fn size(&self) -> usize {
        self.0.len()
    }

    pub fn preprocess<E: EvaluationDomain<Fr>>(
        &self,
        kzg_comm_scheme: &KZGCommitmentScheme,
        domain: &E,
    ) -> PreProcessedTable {
        let coefs = domain.ifft(&self.0);
        let poly = DensePolynomial::from_coefficients_vec(coefs);
        let comm = kzg_comm_scheme.commit(&poly);

        PreProcessedTable {
            poly,
            comm,
            table: self.clone(),
        }
    }

    /// let f = self,
    /// let t = t_table,
    ///
    /// first, we should do is make sure that f is a subset of t,
    /// second, the return table is (f, t) sorted by t.    
    pub fn sort_by(&mut self, t_table: &Self) -> SampleTable {
        let eles = vec![self.0.last().unwrap().to_owned(); t_table.size() - 1 - self.size()];
        self.0.extend(eles);

        let mut sorted = t_table.clone();
        for x in self.0.iter() {
            if !t_table.0.contains(x) {
                panic!("Current prover's table is not a subset of look up table");
            }

            let index = sorted.0.iter().position(|s| *x == *s).unwrap();
            sorted.0.insert(index, *x);
        }

        sorted
    }
}

#[derive(Debug, Clone)]
pub struct LookUpTable {
    t_table: SampleTable,
    f_table: Vec<Fr>,
    domain: Radix2EvaluationDomain<Fr>,
}

impl LookUpTable {
    pub fn new(mut t_table: SampleTable) -> Self {
        if !t_table.size().is_power_of_two() {
            t_table.pad();
        }

        let domain = Radix2EvaluationDomain::new(t_table.size()).unwrap();

        LookUpTable {
            t_table,
            f_table: vec![],
            domain,
        }
    }

    pub fn read_from_u64(&mut self, f: u64) {
        let x = Fr::from(f);
        if !self.f_table.contains(&x) {
            self.f_table.push(x);
        }
    }

    pub fn read_from_scalar(&mut self, f: Fr) {
        if !self.f_table.contains(&f) {
            self.f_table.push(f);
        }
    }

    pub fn prove(&self, kzg_comm_scheme: &KZGCommitmentScheme) -> PlookUpProof {
        for f in self.f_table.iter() {
            if !self.t_table.0.contains(f) {
                panic!("Current prover's table is not a subset of look up table");
            }
        }

        let t_preprocess_table = self.t_table.preprocess(kzg_comm_scheme, &self.domain);
        let f_table = SampleTable::from_scalar(self.f_table.clone());

        prover::prove(f_table, t_preprocess_table, &self.domain, kzg_comm_scheme)
    }
}

pub struct PreProcessedTable {
    pub poly: DensePolynomial<Fr>,
    pub comm: G1Projective,
    pub table: SampleTable,
}
