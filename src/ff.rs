use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, Sub, SubAssign},
    usize,
};

struct ISize {
    value: isize,
}

impl Rem<usize> for ISize {
    type Output = usize;

    fn rem(self, rhs: usize) -> Self::Output {
        if self.value.is_negative() {
            // -a mod b = (b - (a mod b) mod b
            (rhs - (self.value.abs() as usize % rhs)) % rhs
        } else {
            let val = self.value as usize;
            val % rhs
        }
    }
}

pub trait FiniteFieldElement:
    Sized
    + PartialEq
    + Debug
    + Copy
    + Add
    + Sub
    + Mul<Output = Self>
    + Div
    + SubAssign
    + AddAssign
    + MulAssign
    + DivAssign
    + Neg
{
    fn new(element: isize, modulus: usize) -> Self;

    fn element(&self) -> usize;

    fn modulus(&self) -> usize;

    fn zero() -> Self;

    fn one() -> Self;

    fn inverse(&self) -> Self;

    fn pow(&self, n: usize) -> Self;

    fn is_order(&self, order: usize) -> bool;
}

fn multiplicative_inverse(a: isize, b: isize) -> Result<usize, String> {
    /*
     * Computes the multiplicative inverse of a mod b
     * using the "Extended Euclidean Algorithm"
     */

    let modulus = b;

    let (mut m, mut n);
    if a > b {
        m = a;
        n = b;
    } else {
        m = b;
        n = a;
    }
    let mut q = m / n; // quotient
    let mut r = m % n; // remainder
    let mut t_0 = 0;
    let mut t_1 = 1;
    let mut t = t_0 - (t_1 * q);

    while r != 0 {
        (m, n) = (n, r);
        (t_0, t_1) = (t_1, t);
        q = m / n;
        r = m % n;
        t = t_0 - (t_1 * q);
    }

    match n {
        1 => Ok(ISize { value: t_1 } % modulus as usize),
        _ => Err(String::from("Multiplicative inverse does not exist")),
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FFE {
    element: usize,
    modulus: Option<usize>,
}

impl FiniteFieldElement for FFE {
    fn new(value: isize, modulus: usize) -> Self {
        let field_element = ISize { value } % modulus;
        Self {
            element: field_element,
            modulus: Some(modulus),
        }
    }

    fn element(&self) -> usize {
        self.element
    }

    fn modulus(&self) -> usize {
        match self.modulus {
            Some(modulus) => modulus,
            None => 0,
        }
    }

    fn zero() -> Self {
        FFE {
            element: 0,
            modulus: None,
        }
    }

    fn one() -> Self {
        FFE {
            element: 1,
            modulus: None,
        }
    }

    fn inverse(&self) -> Self {
        let inv = multiplicative_inverse(
            self.element.try_into().unwrap(),
            self.modulus.unwrap().try_into().unwrap(),
        )
        .unwrap();
        Self {
            element: inv,
            ..*self
        }
    }

    fn pow(&self, mut n: usize) -> Self {
        let mut current_power = self.to_owned();
        let mut result = Self::one();
        result.modulus = self.modulus;
        while n > 0 {
            if n % 2 == 1 {
                result = result * current_power;
            }
            n = n / 2;
            current_power = current_power * current_power;
        }
        result
    }

    fn is_order(&self, order: usize) -> bool {
        /*
         * Checks the order of an element in the finite field
         * i.e the order of an element is n such that 'a ^ n = e' where
         * a is the element and e is the identity element; in this
         * case 1.
         *
         * The naive approach is the compute the element exponent the order(element ^ order)
         * checking that its equal to the identity and iterate through all the values [1, order - 1]
         * to check that check no other produces the same result as the above(i.e element ^ i = identity)
         *
         * We can perform a simple trick on the second part to make this a bit faster:
         *
         * - Instead of performing exponentiation on every iteration, we multiply by the previous multiple(i.e memoization)
         */

        let mut identity = Self::one();
        identity.modulus = self.modulus;
        let exp = self.pow(order);
        let res = if identity == exp {
            let mut res_inner = true;
            let mut mul = *self;
            for _ in 2..order {
                mul *= *self;
                if mul == identity {
                    res_inner = false;
                    break;
                }
            }
            res_inner
        } else {
            false
        };
        return res;
    }
}

enum Ops {
    ADD,
    MUL,
    SUB,
}

fn perform_mod_operation(op: Ops, a: usize, b: usize, n: usize) -> usize {
    match op {
        Ops::ADD => (a + b) % n,
        Ops::MUL => (a * b) % n,
        Ops::SUB => {
            // TODO: investigate safety of conversion
            let (sub, _) = a.overflowing_sub(b);
            let res = ISize {
                value: sub as isize,
            } % n;
            res
        }
    }
}

impl Add for FFE {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        match (self.modulus, rhs.modulus) {
            (None, None) => Self {
                element: self.element + rhs.element,
                ..self
            },
            (_, _) => {
                let modulus: usize;
                if self.modulus.is_some() {
                    modulus = self.modulus.unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::ADD,
                            self.element,
                            rhs.element,
                            modulus,
                        ),
                        ..self
                    }
                } else {
                    modulus = rhs.modulus.unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::ADD,
                            self.element,
                            rhs.element,
                            modulus,
                        ),
                        ..self
                    }
                }
            }
        }
    }
}

impl AddAssign for FFE {
    fn add_assign(&mut self, rhs: Self) {
        let addition = *self + rhs;
        *self = addition;
    }
}

impl Mul for FFE {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        match (self.modulus, rhs.modulus) {
            (None, None) => Self {
                element: self.element * rhs.element,
                ..self
            },
            (_, _) => {
                let modulus: usize;
                if self.modulus.is_some() {
                    modulus = self.modulus.unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::MUL,
                            self.element,
                            rhs.element,
                            modulus,
                        ),
                        ..self
                    }
                } else {
                    modulus = rhs.modulus.unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::MUL,
                            self.element,
                            rhs.element,
                            modulus,
                        ),
                        ..self
                    }
                }
            }
        }
    }
}

impl MulAssign for FFE {
    fn mul_assign(&mut self, rhs: Self) {
        let mul = *self * rhs;
        *self = mul;
    }
}

impl Sub for FFE {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        match (self.modulus, rhs.modulus) {
            (None, None) => Self {
                element: self.element - rhs.element,
                ..self
            },
            (_, _) => {
                let modulus: usize;
                if self.modulus.is_some() {
                    modulus = self.modulus.unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::SUB,
                            self.element,
                            rhs.element,
                            modulus,
                        ),
                        ..self
                    }
                } else {
                    modulus = rhs.modulus.unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::SUB,
                            self.element,
                            rhs.element,
                            modulus,
                        ),
                        ..self
                    }
                }
            }
        }
    }
}

impl SubAssign for FFE {
    fn sub_assign(&mut self, rhs: Self) {
        let sub = *self - rhs;
        *self = sub;
    }
}

impl Div for FFE {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let inv = rhs.inverse();
        let res = self * inv;
        res
    }
}

impl DivAssign for FFE {
    fn div_assign(&mut self, rhs: Self) {
        let div = *self / rhs;
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

    const MODULUS: usize = 3221225473;

    #[test]
    fn mi() {
        let mi_1 = multiplicative_inverse(3, 5);
        assert_eq!(mi_1, Ok(2));
        let mi_2 = multiplicative_inverse(11, 26);
        assert_eq!(mi_2, Ok(19));
        let mi_2 = multiplicative_inverse(10, 5);
        assert!(mi_2.is_err());
    }

    #[test]
    fn assignment() {
        let ffe_1 = FFE::new(-56, MODULUS);
        let ffe_2 = FFE::new(9704, MODULUS);
        let ffe_3 = FFE::new(3221225477, MODULUS);
        assert_eq!(ffe_1, FFE::new(3221225417, MODULUS));
        assert_eq!(ffe_2, FFE::new(9704, MODULUS));
        assert_eq!(ffe_3, FFE::new(4, MODULUS));
    }

    #[test]
    fn add() {
        let ffe_1 = FFE::new(56, MODULUS);
        let ffe_2 = FFE::new(8902, MODULUS);
        let new_ff = ffe_1 + ffe_2;
        assert_eq!(new_ff, FFE::new(8958, MODULUS));
    }

    #[test]
    fn add_assign() {
        let mut ffe_1 = FFE::new(56, MODULUS);
        let ffe_2 = FFE::new(8902, MODULUS);
        ffe_1 += ffe_2;
        assert_eq!(ffe_1, FFE::new(8958, MODULUS));
    }

    // #[test]
    // fn mul() {
    //     let ffe_1 = FFE::new(1912323, MODULUS);
    //     let ffe_2 = FFE::new(111091, MODULUS);
    //     let new_ff = ffe_1 * ffe_2;
    //     assert_eq!(
    //         new_ff.unwrap(),
    //         FFE {
    //             element: 3062218648,
    //             modulus: MODULUS
    //         }
    //     );
    // }

    // #[test]
    // fn mul_assign() {
    //     let mut ffe_1 = FFE::new(1912323, MODULUS);
    //     let ffe_2 = FFE::new(111091, MODULUS);
    //     ffe_1 *= ffe_2;
    //     assert_eq!(
    //         ffe_1,
    //         FFE {
    //             element: 3062218648,
    //             modulus: MODULUS
    //         }
    //     );
    // }

    // #[test]
    // fn sub() {
    //     let ffe_1 = FFE::new(892, MODULUS);
    //     let ffe_2 = FFE::new(7, MODULUS);
    //     let new_ff = ffe_1 - ffe_2;
    //     assert_eq!(
    //         new_ff.unwrap(),
    //         FFE {
    //             element: 885,
    //             modulus: MODULUS
    //         }
    //     );

    //     let ffe_3 = FFE::new(2, MODULUS);
    //     let ffe_4 = FFE::new(11, MODULUS);
    //     let new_ff = ffe_3 - ffe_4;
    //     assert_eq!(
    //         new_ff.unwrap(),
    //         FFE {
    //             element: 3221225464,
    //             modulus: MODULUS
    //         }
    //     );
    // }

    // #[test]
    // fn sub_assign() {
    //     let mut ffe_1 = FFE::new(2, MODULUS);
    //     let ffe_2 = FFE::new(11, MODULUS);
    //     ffe_1 -= ffe_2;
    //     assert_eq!(
    //         ffe_1,
    //         FFE {
    //             element: 3221225464,
    //             modulus: MODULUS
    //         }
    //     );
    // }

    // #[test]
    // fn div() {
    //     let ffe_1 = FFE::new(892, MODULUS);
    //     let ffe_2 = FFE::new(7, MODULUS);
    //     let new_ff = ffe_1 / ffe_2;
    //     assert_eq!(
    //         new_ff.unwrap(),
    //         FFE {
    //             element: 460175195,
    //             modulus: MODULUS
    //         }
    //     );

    //     let ffe_3 = FFE::new(2, MODULUS);
    //     let ffe_4 = FFE::new(11, MODULUS);
    //     let new_ff = ffe_3 / ffe_4;
    //     assert_eq!(
    //         new_ff.unwrap(),
    //         FFE {
    //             element: 1464193397,
    //             modulus: MODULUS
    //         }
    //     );
    // }

    // #[test]
    // fn div_assign() {
    //     let mut ffe_1 = FFE::new(892, MODULUS);
    //     let ffe_2 = FFE::new(7, MODULUS);
    //     ffe_1 /= ffe_2;
    //     assert_eq!(
    //         ffe_1,
    //         FFE {
    //             element: 460175195,
    //             modulus: MODULUS
    //         }
    //     );

    //     let mut ffe_3 = FFE::new(2, MODULUS);
    //     let ffe_4 = FFE::new(11, MODULUS);
    //     ffe_3 /= ffe_4;
    //     assert_eq!(
    //         ffe_3,
    //         FFE {
    //             element: 1464193397,
    //             modulus: MODULUS
    //         }
    //     );
    // }

    // #[test]
    // fn pow() {
    //     let ffe_1 = FFE::new(76, MODULUS);
    //     let new_ff = ffe_1.pow(2);
    //     assert_eq!(
    //         new_ff,
    //         FFE {
    //             element: 5776,
    //             modulus: MODULUS
    //         }
    //     );

    //     let ffe_2 = FFE::new(700, MODULUS);
    //     let new_ff = ffe_2.pow(90);
    //     assert_eq!(
    //         new_ff,
    //         FFE {
    //             element: 1516783203,
    //             modulus: MODULUS
    //         }
    //     );
    // }

    // #[test]
    // fn is_order() {
    //     #[derive(Debug, Clone, Copy, PartialEq)]
    //     struct SampleFF1 {}

    //     impl FF for SampleFF1 {
    //         type FieldType = usize;
    //         const MODULUS: usize = 71;
    //     }

    //     let ffe_1 = SampleFFE::<SampleFF1>::new(13);
    //     let order: usize = 70;
    //     assert!(ffe_1.is_order(order));
    // }
}
