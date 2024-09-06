use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, Sub, SubAssign},
    u128,
};

struct ISize {
    value: i128,
}

impl Rem<u128> for ISize {
    type Output = u128;

    fn rem(self, rhs: u128) -> Self::Output {
        if self.value.is_negative() {
            // -a mod b = (b - (a mod b) mod b
            (rhs - (self.value.abs() as u128 % rhs)) % rhs
        } else {
            let val = self.value as u128;
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
    fn new(element: i128, modulus: u128) -> Self;

    fn element(&self) -> u128;

    fn modulus(&self) -> u128;

    fn zero() -> Self;

    fn one() -> Self;

    fn inverse(&self) -> Self;

    fn pow(&self, n: u128) -> Self;

    fn is_order(&self, order: u128) -> bool;
}

fn multiplicative_inverse(a: i128, b: i128) -> Result<u128, String> {
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
        1 => Ok(ISize { value: t_1 } % modulus as u128),
        _ => Err(String::from("Multiplicative inverse does not exist")),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct FFE {
    element: u128,
    modulus: Option<u128>,
}

impl FiniteFieldElement for FFE {
    fn new(value: i128, modulus: u128) -> Self {
        let field_element = ISize { value } % modulus;
        Self {
            element: field_element,
            modulus: Some(modulus),
        }
    }

    fn element(&self) -> u128 {
        self.element
    }

    fn modulus(&self) -> u128 {
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

    fn pow(&self, mut n: u128) -> Self {
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

    fn is_order(&self, order: u128) -> bool {
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

fn perform_mod_operation(op: Ops, a: u128, b: u128, n: u128) -> u128 {
    match op {
        Ops::ADD => (a + b) % n,
        Ops::MUL => (a * b) % n,
        Ops::SUB => {
            // TODO: investigate safety of conversion
            let (sub, _) = a.overflowing_sub(b);
            let res = ISize { value: sub as i128 } % n;
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
                let modulus: u128;
                if self.modulus.is_some() {
                    modulus = self.modulus.unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::ADD,
                            self.element,
                            rhs.element,
                            modulus,
                        ),
                        modulus: Some(modulus),
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
                        modulus: Some(modulus),
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
                let modulus: u128;
                if self.modulus.is_some() {
                    modulus = self.modulus.unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::MUL,
                            self.element,
                            rhs.element,
                            modulus,
                        ),
                        modulus: Some(modulus),
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
                        modulus: Some(modulus),
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
            (None, None) => {
                println!("{}", self.element);
                println!("{}", rhs.element);
                Self {
                    element: self.element - rhs.element,
                    ..self
                }
            }
            (_, _) => {
                let modulus: u128;
                if self.modulus.is_some() {
                    modulus = self.modulus.unwrap();
                    Self {
                        element: perform_mod_operation(
                            Ops::SUB,
                            self.element,
                            rhs.element,
                            modulus,
                        ),
                        modulus: Some(modulus),
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
                        modulus: Some(modulus),
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
        match (self.modulus, rhs.modulus) {
            (None, None) => Self {
                element: 0,
                modulus: None,
            },
            (_, _) => {
                if rhs.modulus.is_none() {
                    let rhs = Self {
                        element: rhs.element,
                        ..self
                    };
                    let inv = rhs.inverse();
                    let res = self * inv;
                    return res;
                } else {
                    let inv = rhs.inverse();
                    let res = self * inv;
                    return res;
                }
            }
        }
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

    const MODULUS: u128 = 3221225473;

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
    fn add() {
        let ffe_1 = FFE::new(56, MODULUS);
        let ffe_2 = FFE::new(8902, MODULUS);
        let new_ff = ffe_1 + ffe_2;
        assert_eq!(new_ff, FFE::new(8958, MODULUS));

        let new_ff = ffe_1 + FFE::zero();
        assert_eq!(new_ff, ffe_1);

        let new_ff = FFE::zero() + ffe_1;
        assert_eq!(new_ff, ffe_1);

        let new_ff = ffe_2 + FFE::one();
        assert_eq!(new_ff, FFE::new(8902 + 1, MODULUS));

        let new_ff = FFE::one() + ffe_2;
        assert_eq!(new_ff, FFE::new(8902 + 1, MODULUS));

        let mut ffe_3 = FFE::new(6579, MODULUS);
        let sum = ffe_3 + ffe_1;
        ffe_3 += ffe_1;
        assert_eq!(sum, ffe_3);
    }

    #[test]
    fn mul() {
        let ffe_1 = FFE::new(1912323, MODULUS);
        let ffe_2 = FFE::new(111091, MODULUS);
        let new_ff = ffe_1 * ffe_2;
        assert_eq!(new_ff, FFE::new(3062218648, MODULUS));

        let new_ff = ffe_1 * FFE::zero();
        assert_eq!(new_ff, FFE::new(0, MODULUS));

        let new_ff = FFE::zero() * ffe_1;
        assert_eq!(new_ff, FFE::new(0, MODULUS));

        let new_ff = ffe_1 * FFE::one();
        assert_eq!(new_ff, ffe_1);

        let new_ff = FFE::one() * ffe_1;
        assert_eq!(new_ff, ffe_1);

        let mut ffe_3 = FFE::new(59079, MODULUS);
        ffe_3 *= ffe_1;
        assert_eq!(ffe_3, FFE::new(1912323 * 59079, MODULUS));
    }

    #[test]
    fn sub() {
        let ffe_1 = FFE::new(892, MODULUS);
        let ffe_2 = FFE::new(7, MODULUS);
        let new_ff = ffe_1 - ffe_2;
        assert_eq!(new_ff, FFE::new(885, MODULUS));

        let ffe_3 = FFE::new(2, MODULUS);
        let ffe_4 = FFE::new(11, MODULUS);
        let new_ff = ffe_3 - ffe_4;
        assert_eq!(new_ff, FFE::new(3221225464, MODULUS));

        let new_ff = ffe_1 - FFE::zero();
        assert_eq!(new_ff, ffe_1);

        let new_ff = FFE::zero() - ffe_1;
        assert_eq!(new_ff, FFE::new(-892, MODULUS));

        let new_ff = ffe_1 - FFE::one();
        assert_eq!(new_ff, FFE::new(892 - 1, MODULUS));

        let new_ff = FFE::one() - ffe_1;
        assert_eq!(new_ff, FFE::new(1 - 892, MODULUS));

        let mut ffe_5 = FFE::new(97917, MODULUS);
        ffe_5 -= ffe_1;
        assert_eq!(ffe_5, FFE::new(97917 - 892, MODULUS));
    }

    #[test]
    fn div() {
        let ffe_1 = FFE::new(892, MODULUS);
        let ffe_2 = FFE::new(7, MODULUS);
        let new_ff = ffe_1 / ffe_2;
        assert_eq!(new_ff, FFE::new(460175195, MODULUS));

        let ffe_3 = FFE::new(2, MODULUS);
        let ffe_4 = FFE::new(11, MODULUS);
        let new_ff = ffe_3 / ffe_4;
        assert_eq!(new_ff, FFE::new(1464193397, MODULUS));

        let new_ff = FFE::zero() / ffe_1;
        assert_eq!(new_ff, FFE::new(0, MODULUS));

        let new_ff = ffe_1 / FFE::one();
        assert_eq!(new_ff, FFE::new(892, MODULUS));

        let new_ff = FFE::one() / ffe_1;
        assert_eq!(new_ff, FFE::new(2350916797, MODULUS));

        let mut ffe_5 = FFE::new(892, MODULUS);
        ffe_5 /= ffe_2;
        assert_eq!(ffe_5, FFE::new(460175195, MODULUS));
    }

    #[test]
    fn pow() {
        let ffe_1 = FFE::new(76, MODULUS);
        let new_ff = ffe_1.pow(2);
        assert_eq!(new_ff, FFE::new(5776, MODULUS));

        let ffe_2 = FFE::new(700, MODULUS);
        let new_ff = ffe_2.pow(90);
        assert_eq!(new_ff, FFE::new(1516783203, MODULUS));
    }

    #[test]
    fn is_order() {
        let ffe_1 = FFE::new(13, 71);
        let order: u128 = 70;
        assert!(ffe_1.is_order(order));
    }
}
