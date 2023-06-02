use ark_ff::Field;
use linear_combination::{LinearCombination, Variable};

pub mod constraint_system;
pub mod linear_combination;
pub mod matrices;

pub trait Circuit<F: Field> {
    fn synthesize<R: R1CS<F>>(&self, cs: &mut R);
}

pub trait R1CS<F: Field> {
    fn new_witness(&mut self, value: F) -> Variable<F>;

    fn new_instance(&mut self, value: F) -> Variable<F>;

    fn constrain(
        &mut self,
        a: LinearCombination<F>,
        b: LinearCombination<F>,
        c: LinearCombination<F>,
    );
}
