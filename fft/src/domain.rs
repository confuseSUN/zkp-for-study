use ark_ff::PrimeField;

pub fn get_nth_root_of_unity<F: PrimeField>(size: usize) -> Option<F> {
    let n = size.next_power_of_two();
    let log_n = n.trailing_zeros();

    if n != size || log_n > F::TWO_ADICITY {
        return None;
    }

    let mut omega = F::TWO_ADIC_ROOT_OF_UNITY;
    for _ in log_n..F::TWO_ADICITY {
        omega.square_in_place();
    }

    Some(omega)
}

pub fn evaluation_domain<F: PrimeField>(n: usize) -> Vec<F> {
    assert!(n.is_power_of_two());

    let root: F = get_nth_root_of_unity(n).unwrap();

    let mut elements = Vec::with_capacity(n >> 1);
    elements.extend((0..n >> 1).scan(F::ONE, |state, _| {
        let res = state.clone();
        state.mul_assign(&root);
        Some(res)
    }));

    elements
}
