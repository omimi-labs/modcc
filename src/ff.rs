use num_bigint::BigInt;
use num_traits::One;
use num_traits::Zero;
use std::{
    fmt::{Debug, Display},
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    u128,
};

pub trait FiniteFieldElement:
    Sized
    + Debug
    + Add
    + Sub
    + Mul<Output = Self>
    + Div
    + SubAssign
    + AddAssign
    + MulAssign
    + DivAssign
    + Neg
    + Default
    + PartialEq
    + Eq
    + PartialOrd
    + Ord
{
    fn new(element: &BigInt, modulus: &BigInt) -> Self;

    fn element(&self) -> BigInt;

    fn modulus(&self) -> Option<BigInt>;

    fn zero() -> Self;

    fn one() -> Self;

    fn is_zero(&self) -> bool;

    fn is_one(&self) -> bool;

    fn inverse(&self) -> Option<Self>;

    fn pow(&self, n: u128) -> Self;
}

#[derive(Clone, Default, PartialOrd, Ord)]
pub struct FFE {
    element: BigInt,
    modulus: Option<BigInt>,
}

impl Debug for FFE {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.element)?;
        Ok(())
    }
}

impl PartialEq for FFE {
    fn eq(&self, other: &Self) -> bool {
        if self.is_zero() && other.is_zero() || self.is_one() && other.is_one() {
            true
        } else {
            self.element == other.element && self.modulus == other.modulus
        }
    }
}

impl Eq for FFE {}

impl FiniteFieldElement for FFE {
    fn new(value: &BigInt, modulus: &BigInt) -> Self {
        let field_element = ((value % modulus) + modulus) % modulus;
        Self {
            element: field_element,
            modulus: Some(modulus.clone()),
        }
    }

    fn element(&self) -> BigInt {
        self.element.clone()
    }

    fn modulus(&self) -> Option<BigInt> {
        self.modulus.clone()
    }

    fn zero() -> Self {
        FFE {
            element: BigInt::zero(),
            modulus: None,
        }
    }

    fn one() -> Self {
        FFE {
            element: BigInt::one(),
            modulus: None,
        }
    }

    fn is_zero(&self) -> bool {
        self.element.is_zero()
    }

    fn is_one(&self) -> bool {
        self.element.is_one()
    }

    fn inverse(&self) -> Option<Self> {
        let modulus = self.modulus.clone().unwrap();
        let inv = self.element.modinv(&modulus);
        match inv {
            Some(inverse) => Some(Self {
                element: inverse,
                ..self.clone()
            }),
            None => None,
        }
    }

    fn pow(&self, mut n: u128) -> Self {
        let mut current_power = self.to_owned();
        let mut result = Self::one();
        result.modulus = self.modulus.clone();
        while n > 0 {
            if n % 2 == 1 {
                result = result * current_power.clone();
            }
            n = n / 2;
            current_power *= current_power.clone();
        }
        result
    }
}

enum Ops {
    ADD,
    MUL,
    SUB,
}

fn perform_mod_operation(op: Ops, a: &BigInt, b: &BigInt, n: &BigInt) -> BigInt {
    match op {
        Ops::ADD => ((a + b) + n) % n,
        Ops::MUL => ((a * b) + n) % n,
        Ops::SUB => ((a - b) + n) % n,
    }
}

impl Add for FFE {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        match (&self.modulus, &rhs.modulus) {
            (None, None) => Self {
                element: &self.element + &rhs.element,
                ..self.clone()
            },
            (_, _) => {
                let modulus: BigInt;
                if self.modulus.is_some() {
                    modulus = self.modulus.unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::ADD,
                            &self.element,
                            &rhs.element,
                            &modulus,
                        ),
                        modulus: Some(modulus.clone()),
                    }
                } else {
                    modulus = rhs.modulus.clone().unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::ADD,
                            &self.element,
                            &rhs.element,
                            &modulus,
                        ),
                        modulus: Some(modulus.clone()),
                    }
                }
            }
        }
    }
}

impl AddAssign for FFE {
    fn add_assign(&mut self, rhs: Self) {
        let addition = self.clone() + rhs;
        *self = addition;
    }
}

impl Mul for FFE {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        match (&self.modulus, &rhs.modulus) {
            (None, None) => Self {
                element: self.element.clone() * rhs.element,
                ..self.clone()
            },
            (_, _) => {
                let modulus: BigInt;
                if self.modulus.is_some() {
                    modulus = self.clone().modulus.clone().unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::MUL,
                            &self.element,
                            &rhs.element,
                            &modulus,
                        ),
                        modulus: Some(modulus.clone()),
                    }
                } else {
                    modulus = rhs.modulus.unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::MUL,
                            &self.element,
                            &rhs.element,
                            &modulus,
                        ),
                        modulus: Some(modulus),
                    }
                }
            }
        }
    }
}

impl MulAssign for FFE {
    fn mul_assign(&mut self, rhs: Self) {
        let mul = self.clone() * rhs;
        *self = mul;
    }
}

impl Sub for FFE {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        match (&self.modulus, &rhs.modulus) {
            (None, None) => Self {
                element: &self.element - &rhs.element,
                ..self
            },
            (_, _) => {
                let modulus: BigInt;
                if self.modulus.is_some() {
                    modulus = self.modulus.unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::SUB,
                            &self.element,
                            &rhs.element,
                            &modulus,
                        ),
                        modulus: Some(modulus),
                    }
                } else {
                    modulus = rhs.modulus.unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::SUB,
                            &self.element,
                            &rhs.element,
                            &modulus,
                        ),
                        modulus: Some(modulus),
                    }
                }
            }
        }
    }
}

impl SubAssign for FFE {
    fn sub_assign(&mut self, rhs: Self) {
        let sub = self.clone() - rhs;
        *self = sub;
    }
}

impl Div for FFE {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        match (&self.modulus, &rhs.modulus) {
            (None, None) => Self {
                element: BigInt::zero(),
                modulus: None,
            },
            (_, _) => {
                if rhs.modulus.is_none() {
                    let rhs = Self {
                        element: rhs.element,
                        ..self.clone()
                    };
                    let inv = rhs.inverse().unwrap();
                    let res = self * inv;
                    return res;
                } else {
                    let inv = rhs.inverse().unwrap();
                    let res = self * inv;
                    return res;
                }
            }
        }
    }
}

impl DivAssign for FFE {
    fn div_assign(&mut self, rhs: Self) {
        let div = self.clone() / rhs;
        *self = div;
    }
}

impl Neg for FFE {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let neg = Self::zero()
            - Self {
                element: self.element,
                ..self
            };
        neg
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MODULUS: u128 = 3221225473;

    #[test]
    fn add() {
        let modulus = BigInt::from(MODULUS);

        let ffe_1 = FFE::new(&BigInt::from(56), &modulus);
        let ffe_2 = FFE::new(&BigInt::from(8902), &modulus);
        let new_ff = ffe_1.clone() + ffe_2.clone();
        assert_eq!(new_ff, FFE::new(&BigInt::from(8958), &modulus));

        let new_ff = ffe_1.clone() + FFE::zero();
        assert_eq!(new_ff, ffe_1);

        let new_ff = FFE::zero() + ffe_1.clone();
        assert_eq!(new_ff, ffe_1);

        let new_ff = ffe_2.clone() + FFE::one();
        assert_eq!(new_ff, FFE::new(&BigInt::from(8902 + 1), &modulus));

        let new_ff = FFE::one() + ffe_2;
        assert_eq!(new_ff, FFE::new(&BigInt::from(8902 + 1), &modulus));

        let mut ffe_3 = FFE::new(&BigInt::from(6579), &modulus);
        let sum = ffe_3.clone() + ffe_1.clone();
        ffe_3 += ffe_1;
        assert_eq!(sum, ffe_3);
    }

    #[test]
    fn mul() {
        let modulus = BigInt::from(MODULUS);

        let ffe_1 = FFE::new(&BigInt::from(1912323_u32), &modulus);
        let ffe_2 = FFE::new(&BigInt::from(111091_u32), &modulus);
        let new_ff = ffe_1.clone() * ffe_2;
        assert_eq!(new_ff, FFE::new(&BigInt::from(3062218648_u32), &modulus));

        let new_ff = ffe_1.clone() * FFE::zero();
        assert_eq!(new_ff, FFE::new(&BigInt::zero(), &modulus));

        let new_ff = FFE::zero() * ffe_1.clone();
        assert_eq!(new_ff, FFE::new(&BigInt::zero(), &modulus));

        let new_ff = ffe_1.clone() * FFE::one();
        assert_eq!(new_ff, ffe_1);

        let new_ff = FFE::one() * ffe_1.clone();
        assert_eq!(new_ff, ffe_1);

        let mut ffe_3 = FFE::new(&BigInt::from(59079), &modulus);
        ffe_3 *= ffe_1;
        assert_eq!(
            ffe_3,
            FFE::new(&BigInt::from(1912323_u64 * 59079_u64), &modulus)
        );
    }

    #[test]
    fn sub() {
        let modulus = BigInt::from(MODULUS);

        let ffe_1 = FFE::new(&BigInt::from(892), &modulus);
        let ffe_2 = FFE::new(&BigInt::from(7), &modulus);

        let new_ff = ffe_1.clone() - ffe_2;
        assert_eq!(new_ff, FFE::new(&BigInt::from(885), &modulus));

        let ffe_3 = FFE::new(&BigInt::from(2), &modulus);
        let ffe_4 = FFE::new(&BigInt::from(11), &modulus);
        let new_ff = ffe_3 - ffe_4;
        assert_eq!(new_ff, FFE::new(&BigInt::from(3221225464_u32), &modulus));

        let new_ff = ffe_1.clone() - FFE::zero();
        assert_eq!(new_ff, ffe_1.clone());

        let new_ff = FFE::zero() - ffe_1.clone();
        assert_eq!(new_ff, FFE::new(&BigInt::from(-892), &modulus));

        let new_ff = ffe_1.clone() - FFE::one();
        assert_eq!(new_ff, FFE::new(&BigInt::from(892 - 1), &modulus));

        let new_ff = FFE::one() - ffe_1.clone();
        assert_eq!(new_ff, FFE::new(&BigInt::from(1 - 892), &modulus));

        let mut ffe_5 = FFE::new(&BigInt::from(97917), &modulus);
        ffe_5 -= ffe_1;
        assert_eq!(ffe_5, FFE::new(&BigInt::from(97917 - 892), &modulus));
    }

    #[test]
    fn div() {
        let modulus = BigInt::from(MODULUS);

        let ffe_1 = FFE::new(&BigInt::from(892), &modulus);
        let ffe_2 = FFE::new(&BigInt::from(7), &modulus);
        let new_ff = ffe_1.clone() / ffe_2.clone();
        assert_eq!(new_ff, FFE::new(&BigInt::from(460175195), &modulus));

        let ffe_3 = FFE::new(&BigInt::from(2), &modulus);
        let ffe_4 = FFE::new(&BigInt::from(11), &modulus);
        let new_ff = ffe_3 / ffe_4;
        assert_eq!(new_ff, FFE::new(&BigInt::from(1464193397), &modulus));

        let new_ff = FFE::zero() / ffe_1.clone();
        assert_eq!(new_ff, FFE::new(&BigInt::zero(), &modulus));

        let new_ff = ffe_1.clone() / FFE::one();
        assert_eq!(new_ff, FFE::new(&BigInt::from(892), &modulus));

        let new_ff = FFE::one() / ffe_1.clone();
        assert_eq!(new_ff, FFE::new(&BigInt::from(2350916797_u32), &modulus));

        let mut ffe_5 = FFE::new(&BigInt::from(892), &modulus);
        ffe_5 /= ffe_2;
        assert_eq!(ffe_5, FFE::new(&BigInt::from(460175195), &modulus));
    }

    #[test]
    fn no_modulus() {
        let zero = FFE::zero();
        let one = FFE::one();

        // 1 + 1
        let ffe = one.clone() + one.clone();
        assert_eq!(
            ffe,
            FFE {
                element: BigInt::from(2),
                modulus: None
            }
        );
        // 0 + 0
        let ffe = zero.clone() + zero.clone();
        assert_eq!(
            ffe,
            FFE {
                element: BigInt::zero(),
                modulus: None
            }
        );
        // 1 - 0
        let ffe = one.clone() - zero.clone();
        assert_eq!(
            ffe,
            FFE {
                element: BigInt::one(),
                modulus: None
            }
        );
        // 0 - 1
        let ffe = zero - one;
        assert_eq!(
            ffe,
            FFE {
                element: BigInt::from(-1),
                modulus: None
            }
        )
    }

    #[test]
    fn pow() {
        let modulus = BigInt::from(MODULUS);

        let ffe_1 = FFE::new(&BigInt::from(76), &modulus);
        let new_ff = ffe_1.pow(2);
        assert_eq!(new_ff, FFE::new(&BigInt::from(5776), &modulus));

        let ffe_2 = FFE::new(&BigInt::from(700), &modulus);
        let new_ff = ffe_2.pow(90);
        assert_eq!(new_ff, FFE::new(&BigInt::from(1516783203), &modulus));
    }
}
