pub fn butterfly<T>(input: &mut [T]) {
    assert!(input.len().is_power_of_two());

    let n = input.len();
    let log_n = n.trailing_zeros();

    let reverse_bits = |i: usize| -> usize {
        if n == 1 {
            i
        } else {
            i.reverse_bits() >> (usize::BITS - log_n)
        }
    };

    for i in 0..n {
        let rev = reverse_bits(i);
        if rev > i {
            input.swap(i, rev)
        }
    }
}

#[cfg(test)]
mod test {
    use super::butterfly;

    #[test]
    fn test_butterfly() {
        let x = vec![0, 1, 2, 3, 4, 5, 6, 7];
        let mut y = x.clone();

        butterfly(&mut y);
        assert_eq!(y, vec![0, 4, 2, 6, 1, 5, 3, 7]);

        butterfly(&mut y);
        assert_eq!(x, y);
    }
}
