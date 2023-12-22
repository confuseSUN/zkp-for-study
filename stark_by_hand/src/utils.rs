use std::ops::AddAssign;

use ark_ff::PrimeField;

use crate::column::Column;

// pub fn vec_mul<F:PrimeField>(a:&[F],b:&[F]) -> Vec<F> {
//     assert_eq!(a.len(),b.len());

//      let mut c = Vec::new();
//     for (x,y)  in a.iter().zip(b.iter()) {

//     }
// }


pub fn mix<F:PrimeField>(columns : &[Column<F>],alpha :&F)  -> Column<F>{
    let mut res = columns[0].clone();

    let mut multiplier = alpha.clone();
    
    for c in columns.iter().skip(1) {
        res.add_assign(&c.scale(&multiplier));
        multiplier.mul_assign(alpha);    
    }

    res

}
