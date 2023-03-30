use ark_ff::fields::{Fp256, Fp64, MontBackend, MontConfig};

#[derive(MontConfig)]
#[modulus = "101"]
#[generator = "2"]
pub struct FrConfig;
pub type F101 = Fp64<MontBackend<FrConfig, 1>>;

#[derive(MontConfig)]
#[modulus = "21888242871839275222246405745257275088548364400416034343698204186575808495617"]
#[generator = "7"]
pub struct BN254FrConfig;
pub type BN254Fr = Fp256<MontBackend<BN254FrConfig, 4>>;
