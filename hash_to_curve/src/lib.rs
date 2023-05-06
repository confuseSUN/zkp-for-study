use ark_ec::hashing::curve_maps::swu::parity;
use ark_ec::short_weierstrass::Affine;
use ark_ec::short_weierstrass::SWCurveConfig;
use ark_ec::AffineRepr;
use ark_ff::field_hashers::{DefaultFieldHasher, HashToField};
use ark_ff::Field;
use sha2::digest::DynDigest;
use std::ops::*;

/// We have implemented a general Shallue–van de Woestijne map that is effective for most short Weierstrass equation g(x) = y^2 = x^3 + a * x +b.
///
/// To implement hashing to curve, the following steps are usually required,
/// see <https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-16.html> :
/// step 1 : u = hash_to_field(msg, 2)
/// step 2 : Q0 = map_to_curve(u[0])
/// step 3 : Q1 = map_to_curve(u[1])
/// step 4 : R = Q0 + Q1
/// step 5 : P = clear_cofactor(R)
/// step 6 : return P
pub trait HashToCurve<P: SWCurveConfig> {
    /// Z is a non-zero element meeting the below criteria:
    /// 1. g(Z) != 0
    /// 2. -(3 * Z^2 + 4 * A) / (4 * g(Z)) != 0
    /// 3. -(3 * Z^2 + 4 * A) / (4 * g(Z)) is square
    /// 4. at least one of g(Z) and g(-Z / 2) is square
    const Z: P::BaseField;

    /// Domain separation
    const DST: &'static [u8];

    /// Mapping an arbitrary field element to a point on the elliptic curve,
    /// This step matching step 2 and step 3
    fn map_to_curve(point: &P::BaseField) -> Result<Affine<P>, String> {
        let a: P::BaseField = P::COEFF_A;
        let b: P::BaseField = P::COEFF_B;

        let tv1: P::BaseField = point.mul(point).mul(&Self::c1());
        let tv2: P::BaseField = P::BaseField::ONE.add(&tv1);
        let tv1: P::BaseField = P::BaseField::ONE.sub(&tv1);
        let tv3: P::BaseField = tv1.mul(&tv2);
        let tv3: P::BaseField = tv3.inverse().unwrap();
        let tv4: P::BaseField = point.mul(&tv1).mul(&tv3).mul(&Self::c3());

        let x1: P::BaseField = Self::c2().sub(&tv4);
        let gx1: P::BaseField = x1.mul(&x1).add(&a);
        let gx1: P::BaseField = gx1.mul(&x1).add(&b);
        if gx1.legendre().is_qr() {
            let y: P::BaseField = gx1.sqrt().unwrap();
            let y: P::BaseField = if parity(&y) != parity(point) { -y } else { y };
            let point_on_curve = Affine::<P>::new_unchecked(x1, y);
            assert!(point_on_curve.is_on_curve());
            return Ok(point_on_curve);
        }

        let x2: P::BaseField = Self::c2().add(&tv4);
        let gx2: P::BaseField = x2.mul(&x2).add(&a);
        let gx2: P::BaseField = gx2.mul(&x2).add(&b);
        if gx2.legendre().is_qr() {
            let y: P::BaseField = gx2.sqrt().unwrap();
            let y: P::BaseField = if parity(&y) != parity(point) { -y } else { y };
            let point_on_curve = Affine::<P>::new_unchecked(x2, y);
            assert!(point_on_curve.is_on_curve());
            return Ok(point_on_curve);
        }

        let x3: P::BaseField = tv2.mul(&tv2).mul(&tv3);
        let x3: P::BaseField = x3.mul(&x3).mul(&Self::c4()).add(&Self::Z);
        let gx3: P::BaseField = x3.mul(&x3).add(&a);
        let gx3: P::BaseField = gx3.mul(&x3).add(&b);
        if gx3.legendre().is_qr() {
            let y: P::BaseField = gx3.sqrt().unwrap();
            let y: P::BaseField = if parity(&y) != parity(point) { -y } else { y };
            let point_on_curve = Affine::<P>::new_unchecked(x3, y);
            assert!(point_on_curve.is_on_curve());
            return Ok(point_on_curve);
        }

        Err("faile to hash to curve".to_string())
    }

    /// Mapping an arbitrary message to a field element,
    /// This step matching step 1
    fn hash_to_field<H: Default + DynDigest + Clone>(msg: &[u8], dst: &[u8]) -> Vec<P::BaseField> {
        let field_hasher = <DefaultFieldHasher<H> as HashToField<P::BaseField>>::new(dst);

        let field_elems: Vec<P::BaseField> = field_hasher.hash_to_field(msg, 2);

        field_elems
    }

    /// Perform hashing to curve
    fn hash<H: Default + DynDigest + Clone>(msg: &[u8]) -> Result<Affine<P>, String> {
        let rand_field_elems: Vec<P::BaseField> = Self::hash_to_field::<H>(msg, Self::DST);

        let rand_curve_elem_0 = Self::map_to_curve(&rand_field_elems[0])?;
        let rand_curve_elem_1 = Self::map_to_curve(&rand_field_elems[1])?;

        let rand_curve_elem: Affine<P> = rand_curve_elem_0.add(&rand_curve_elem_1).into();
        let rand_subgroup_elem = rand_curve_elem.clear_cofactor();

        Ok(rand_subgroup_elem)
    }

    /// The constant c1 equals ：
    /// c1 = g(Z)
    fn c1() -> P::BaseField {
        let a = P::COEFF_A;
        let b = P::COEFF_B;

        let c1 = Self::Z.mul(&Self::Z).add(&a);
        let c1 = c1.mul(&Self::Z).add(&b);

        c1
    }

    /// The constant c2 equals ：
    /// c2 = -Z / 2
    fn c2() -> P::BaseField {
        let two = P::BaseField::ONE.add(&P::BaseField::ONE);
        let c2 = Self::Z.div(two);
        c2.neg()
    }

    /// The constant c3 equals ：
    /// c3 = sqrt(-g(Z) * (3 * Z^2 + 4 * A)) and sgn0(c3) MUST equal 0
    fn c3() -> P::BaseField {
        let a = P::COEFF_A;
        let three = P::BaseField::from(3u64);
        let four = P::BaseField::from(4u64);

        let gz = Self::c1();
        let gz_neg = gz.neg();

        let tmp = Self::Z
            .mul(&Self::Z)
            .mul(&three)
            .add(a.mul(&four))
            .mul(&gz_neg);

        let mut c3 = tmp.sqrt().unwrap();

        if parity(&c3) {
            c3.neg_in_place();
        }

        c3
    }

    /// The constant c4 equals ：
    /// c4 = -4 * g(Z) / (3 * Z^2 + 4 * A)
    fn c4() -> P::BaseField {
        let a = P::COEFF_A;
        let three = P::BaseField::from(3u64);
        let four = P::BaseField::from(4u64);

        let gz = Self::c1();
        let gz_neg = gz.neg();

        let c4 = gz_neg
            .mul(&four)
            .div(Self::Z.mul(&Self::Z).mul(&three).add(a.mul(&four)));
        c4
    }
}

#[cfg(test)]
mod tests {
    use crate::HashToCurve;
    use ark_ec::short_weierstrass::SWCurveConfig;
    use ark_ff::Field;
    use ark_ff::MontFp;
    use ark_ff::Zero;
    use ark_secp256k1::Affine;
    use ark_secp256k1::{Config, Fq};
    use sha2::Sha256;
    use std::ops::Neg;
    use std::ops::{Add, Div, Mul};

    pub struct Secp256K1SWMap;

    impl HashToCurve<ark_secp256k1::Config> for Secp256K1SWMap {
        const Z: Fq = MontFp!("1");

        const DST: &'static [u8] = b"secp256k1";
    }

    #[test]
    fn test_hash_to_curve_for_secp256k1() {
        let msg = b"hello, hash to secp256k1 ";
        let point = Secp256K1SWMap::hash::<Sha256>(msg).unwrap();

        let x = MontFp!(
            "64043297498491436634439171231145164731426160760324877170764636047178421698471"
        );
        let y = MontFp!(
            "91198262659971141530092929731991797815001480798760841105513408186364380289421"
        );
        let expect_point = Affine::new_unchecked(x, y);

        assert!(expect_point.is_on_curve());

        assert_eq!(point, expect_point);
    }

    #[test]
    fn cal_z_for_secp256k1() {
        let candidate = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5];

        let a = Config::COEFF_A;
        let b = Config::COEFF_B;
        let two = Fq::from(2);
        let three = Fq::from(3);
        let four = Fq::from(4);

        let four_a = four.mul(&a);

        for x in candidate {
            let z = Fq::from(x);

            // require g(Z) != 0
            let y = z.mul(&z).add(&a);
            let y = y.mul(&z).add(&b);
            if y == Fq::zero() {
                continue;
            }

            // require -(3 * Z^2 + 4 * A) / (4 * g(Z)) != 0  and -(3 * Z^2 + 4 * A) / (4 * g(Z)) is square
            let four_y = four.mul(&y);
            let four_y_inv = four_y.inverse().unwrap();
            let tmp = z.mul(&z).mul(&three).add(&four_a);
            let tmp = tmp.neg();
            let tmp = tmp.mul(&four_y_inv);
            if tmp == Fq::zero() || !tmp.legendre().is_qr() {
                continue;
            }

            // require at least one of g(Z) and g(-Z / 2) is square
            let z_neg = z.neg();
            let z_neg_div_two = z_neg.div(&two);
            let y1 = z_neg_div_two.mul(&z_neg_div_two).add(&a);
            let y1 = y1.mul(&z_neg_div_two).add(&b);

            if !y.legendre().is_qr() && !y1.legendre().is_qr() {
                continue;
            }

            println!("z: {}", x);
        }
    }
}
