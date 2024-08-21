use crate::ff::FiniteFieldElement;
use std::ops::{Add, Mul, Neg, Sub};

pub trait Polynomial<F>: Sized {
    fn is_zero(&self) -> bool;

    fn degree(&self) -> usize;

    fn new(coefficients: Vec<F>) -> Self;

    fn x() -> Self;

    fn zero() -> Self;

    fn one() -> Self;
}

pub trait UnivariatePolynomial<F>: Polynomial<F> {
    fn coefficients(&self) -> Vec<F>;

    fn evaluate(&self, x: F) -> F;

    fn interpolate(y_values: &Vec<F>) -> Self;

    fn interpolate_xy(x_values: &Vec<F>, y_values: &Vec<F>) -> Self;

    fn get_lagrange_polynomial(x_value: F, x_values: &Vec<F>) -> Self;
}

// Univariate Polynomial
#[derive(Debug, PartialEq, Clone)]
pub struct UniPoly<F> {
    // Co-efficient represented from lower degree to higher
    // For example: 2x^2 + x + 1 is represented as [1, 1, 2]
    coefficients: Vec<F>,
}

impl<F: FiniteFieldElement> Polynomial<F> for UniPoly<F> {
    fn new(coefficients: Vec<F>) -> Self {
        UniPoly { coefficients }
    }

    fn x() -> Self {
        Self::new(vec![F::zero(), F::one()])
    }

    fn one() -> Self {
        Self::new(vec![F::one()])
    }

    fn zero() -> Self {
        Self::new(vec![F::zero()])
    }

    fn is_zero(&self) -> bool {
        self.coefficients.is_empty()
            || (self.coefficients.len() == 1 && self.coefficients[0] == F::zero())
    }

    fn degree(&self) -> usize {
        if self.is_zero() {
            0
        } else {
            self.coefficients.len() - 1
        }
    }
}

impl<F: FiniteFieldElement + Neg<Output = F> + Sub<Output = F> + Add<Output = F>>
    UnivariatePolynomial<F> for UniPoly<F>
{
    fn coefficients(&self) -> Vec<F> {
        self.coefficients.clone()
    }

    fn evaluate(&self, var: F) -> F {
        let mut identity = F::zero();
        for (i, x) in self.coefficients.iter().enumerate() {
            let exp = var.pow(i);
            let mul = exp * *x;
            identity += mul
        }
        identity
    }

    fn interpolate(y_values: &Vec<F>) -> Self {
        let mut x_values = vec![];
        for (i, y) in y_values.iter().enumerate() {
            x_values.push(F::new(i.try_into().unwrap(), y.modulus()));
        }
        Self::interpolate_xy(&x_values, y_values)
    }

    fn interpolate_xy(x_values: &Vec<F>, y_values: &Vec<F>) -> Self {
        assert_eq!(x_values.len(), y_values.len());
        let mut resulting_polynomial = Self::zero();
        for (x, y) in x_values.iter().zip(y_values.iter()) {
            let lagrange_polynomial = Self::get_lagrange_polynomial(*x, &x_values);
            println!("{:?}", lagrange_polynomial);
            let y_poly = Self::new(vec![*y]);
            let product = &y_poly * &lagrange_polynomial;
            resulting_polynomial = &resulting_polynomial + &product;
        }
        resulting_polynomial
    }

    fn get_lagrange_polynomial(x_value: F, x_values: &Vec<F>) -> Self {
        /*
         * L_i = \prod_{j=0, j \neq i}^{n} \dfrac{x - x_j}{x_i - x_j}
         *
         *  where:
         *      `i` is x_value
         *      `j` is the index in the loop below
         */
        let mut resulting_polynomial = Self::one();
        for x in x_values.iter() {
            if *x == x_value {
                continue;
            }
            let numerator = Self::new(vec![-(*x), F::one()]);
            let inverse_of_denominator = Self::new(vec![(x_value - *x).inverse()]);
            let product = &numerator * &inverse_of_denominator;
            // TODO: implement mul_assign() for &UniPoly
            resulting_polynomial = &resulting_polynomial * &product;
        }
        resulting_polynomial
    }
}

impl<F: FiniteFieldElement> Mul for &UniPoly<F> {
    type Output = UniPoly<F>;

    fn mul(self, other: Self) -> Self::Output {
        if self.is_zero() || other.is_zero() {
            UniPoly::zero()
        } else {
            let deg_a = self.coefficients.len() - 1;
            let deg_b = other.coefficients.len() - 1;
            let max_coefficients = deg_a + deg_b + 1;
            let zero = F::zero();
            let mut product_coefficients = vec![zero; max_coefficients];
            for i in 0..=self.degree() {
                for j in 0..=other.degree() {
                    let index = i + j;
                    let product = self.coefficients[i] * other.coefficients[j];
                    product_coefficients[index] += product;
                }
            }
            let poly = UniPoly::new(product_coefficients);
            poly
        }
    }
}

impl<F: FiniteFieldElement + Add<Output = F>> Add for &UniPoly<F> {
    type Output = UniPoly<F>;

    fn add(self, rhs: Self) -> Self::Output {
        if self.is_zero() {
            rhs.clone()
        } else if rhs.is_zero() {
            self.clone()
        } else {
            let new_coefficients = add_list(self.coefficients.clone(), rhs.coefficients.clone());
            UniPoly::new(new_coefficients)
        }
    }
}

fn add_list<T: FiniteFieldElement + Add<Output = T>>(a: Vec<T>, b: Vec<T>) -> Vec<T> {
    let mut res: Vec<T> = vec![];
    if a.len() == b.len() {
        for (x, y) in a.iter().zip(b.iter()) {
            res.push(*x + *y);
        }
    } else if a.len() > b.len() {
        let diff = a.len() - b.len();
        let mut b = b;
        for _ in 0..diff {
            b.push(T::zero());
        }
        for (x, y) in a.iter().zip(b.iter()) {
            res.push(*x + *y);
        }
    } else {
        let diff = b.len() - a.len();
        let mut a = a;
        for _ in 0..diff {
            a.push(T::zero());
        }
        for (x, y) in a.iter().zip(b.iter()) {
            res.push(*x + *y);
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ff::FFE;

    const MODULUS: usize = 3221225473;

    #[test]
    fn mul() {
        // Test:
        let x = UniPoly::<FFE>::x();
        let co_effs = vec![FFE::one(), FFE::one()];
        let x_plus_1 = UniPoly::<FFE>::new(co_effs);
        let actual = &x * &x_plus_1;
        let co_effs = vec![FFE::zero(), FFE::one(), FFE::one()];
        let expected = UniPoly::<FFE>::new(co_effs);
        assert_eq!(actual, expected);

        // Tests (x^3 - 3x + 2) * (2x + 5)
        let co_effs_1 = vec![
            FFE::new(2, MODULUS),
            FFE::new(-3, MODULUS),
            FFE::zero(),
            FFE::new(1, MODULUS),
        ];
        let poly_1 = UniPoly::<FFE>::new(co_effs_1);
        let co_effs_2 = vec![FFE::new(5, MODULUS), FFE::new(2, MODULUS)];
        let poly_2 = UniPoly::<FFE>::new(co_effs_2);
        let actual = &poly_1 * &poly_2;
        let exp_co_effs = vec![
            FFE::new(10, MODULUS),
            FFE::new(-11, MODULUS),
            FFE::new(-6, MODULUS),
            FFE::new(5, MODULUS),
            FFE::new(2, MODULUS),
        ];
        let expected = UniPoly::<FFE>::new(exp_co_effs);
        assert_eq!(actual, expected);
    }

    #[test]
    fn add() {
        // Test: (x^2 + x + 5) + (2x^2 + 4x + 2)
        let co_eff_1 = vec![
            FFE::new(5, MODULUS),
            FFE::new(1, MODULUS),
            FFE::new(1, MODULUS),
        ];
        let poly_1 = UniPoly::new(co_eff_1);
        let co_eff_2 = vec![
            FFE::new(2, MODULUS),
            FFE::new(4, MODULUS),
            FFE::new(2, MODULUS),
        ];
        let poly_2 = UniPoly::new(co_eff_2);
        let actual = &poly_1 + &poly_2;
        let exp_co_effs = vec![
            FFE::new(7, MODULUS),
            FFE::new(5, MODULUS),
            FFE::new(3, MODULUS),
        ];
        let expected = UniPoly::<FFE>::new(exp_co_effs);
        assert_eq!(actual, expected);

        // Test: (x^3 - 3x + 2) + (2x + 5)
        let co_eff_3 = vec![FFE::new(5, MODULUS), FFE::new(2, MODULUS)];
        let poly_3 = UniPoly::new(co_eff_3);
        let co_eff_4 = vec![
            FFE::new(2, MODULUS),
            FFE::new(-3, MODULUS),
            FFE::new(0, MODULUS),
            FFE::new(1, MODULUS),
        ];
        let poly_4 = UniPoly::new(co_eff_4);
        let actual = &poly_3 + &poly_4;
        let exp_co_effs = vec![
            FFE::new(7, MODULUS),
            FFE::new(-1, MODULUS),
            FFE::new(0, MODULUS),
            FFE::new(1, MODULUS),
        ];
        let expected = UniPoly::<FFE>::new(exp_co_effs);
        assert_eq!(actual, expected);
    }

    #[test]
    fn interpolate() {
        // Interpolating the values: [3, 1, 2, 4]
        let modulus = 17;
        let co_effs = vec![
            FFE::new(3, modulus),
            FFE::new(1, modulus),
            FFE::new(2, modulus),
            FFE::new(4, modulus),
        ];
        let polynomial: UniPoly<FFE> = UniPoly::<FFE>::interpolate(&co_effs);
        println!("{:?}", polynomial);
    }
}
