use crate::table::{LookUpTable, SampleTable};
use ark_std::test_rng;
use study_kzg::{commitment::KZGCommitmentScheme, srs::SRS};

#[test]
fn test_plookup() {
    let x = vec![1, 5, 7, 20, 21, 24, 56, 100];
    let t = SampleTable::from_u64(x);
    let mut look = LookUpTable::new(t);
    look.read_from_u64(24);
    look.read_from_u64(21);
    look.read_from_u64(56);
    look.read_from_u64(56);
    look.read_from_u64(56);

    let mut rng = test_rng();
    let max_degree = 20;
    let srs = SRS::new(max_degree, &mut rng);
    let kzg_comm_scheme = KZGCommitmentScheme(&srs);
    let proof = look.prove(&kzg_comm_scheme);
    assert!(proof.verify(&kzg_comm_scheme, &mut rng));
}
