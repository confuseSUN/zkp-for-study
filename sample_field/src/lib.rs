use ark_ff::fields::{Fp64, MontBackend, MontConfig};

#[derive(MontConfig)]
#[modulus = "101"]
#[generator = "2"]
pub struct FrConfig;

pub type F101 = Fp64<MontBackend<FrConfig, 1>>;
