use ark_poly::{EvaluationDomain, Radix2EvaluationDomain, UVPolynomial};
use study_kzg::{commitment::KZGCommitmentScheme, Poly, Scalar, G1};

use crate::{prover, verifier::PlookUpProof};

#[derive(Debug, Clone)]
pub struct SampleTable(pub Vec<Scalar>);

impl SampleTable {
    pub fn from_scalar(t: Vec<Scalar>) -> Self {
        Self(t)
    }

    pub fn from_u64(t: Vec<u64>) -> Self {
        let t = t.into_iter().map(|x| Scalar::from(x)).collect();
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

    pub fn preprocess<E: EvaluationDomain<Scalar>>(
        &self,
        kzg_comm_scheme: &KZGCommitmentScheme,
        domain: &E,
    ) -> PreProcessedTable {
        let coefs = domain.ifft(&self.0);
        let poly = Poly::from_coefficients_vec(coefs);
        let comm = kzg_comm_scheme.commit(&poly);

        PreProcessedTable {
            poly,
            comm,
            table: self.clone(),
        }
    }

    pub fn sort_by(&mut self, sorted_table: &Self) -> SampleTable {
        let eles = vec![self.0.last().unwrap().to_owned(); sorted_table.size() - 1 - self.size()];
        self.0.extend(eles);

        let mut sorted = sorted_table.clone();
        for x in self.0.iter() {
            if !sorted_table.0.contains(x) {
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
    f_table: Vec<Scalar>,
    domain: Radix2EvaluationDomain<Scalar>,
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
        self.f_table.push(Scalar::from(f));
    }

    pub fn read_from_scalar(&mut self, scalar: Scalar) {
        self.f_table.push(scalar);
    }

    pub fn prove(&self, kzg_comm_scheme: &KZGCommitmentScheme) -> PlookUpProof {
        let t_preprocess_table = self.t_table.preprocess(kzg_comm_scheme, &self.domain);
        let f_table = SampleTable::from_scalar(self.f_table.clone());

        prover::prove(f_table, t_preprocess_table, &self.domain, kzg_comm_scheme)
    }
}

pub struct PreProcessedTable {
    pub poly: Poly,
    pub comm: G1,
    pub table: SampleTable,
}
