use ark_ff::{Fp64, MontBackend, MontConfig};

#[derive(MontConfig)]
#[modulus = "97"]
#[generator = "5"]
#[small_subgroup_base = "3"]
#[small_subgroup_power = "1"]
pub struct FrConfig;
pub type Fr = Fp64<MontBackend<FrConfig, 1>>;
