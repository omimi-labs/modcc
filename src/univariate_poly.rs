use std::{
    collections,
    ops::{Add, Mul, Neg, Sub},
    u128,
};

use crate::ff::FiniteFieldElement;
use num_bigint::BigInt;
use serde::Serialize;

pub trait Polynomial<F>: Sized {
    fn is_zero(&self) -> bool;

    fn degree(&self) -> usize;

    fn new(coefficients: Vec<F>) -> Self;

    fn x() -> Self;

    fn zero() -> Self;

    fn one() -> Self;
}

pub trait UnivariatePolynomial<F: Default>: Polynomial<F> {
    fn coefficients(&self) -> Vec<F>;

    fn evaluate(&self, x: &F) -> F;

    fn interpolate(y_values: &Vec<F>) -> (Self, LagrangeInterpolationSteps);

    fn interpolate_xy(x_values: &Vec<F>, y_values: &Vec<F>) -> (Self, LagrangeInterpolationSteps);

    fn get_lagrange_polynomial(x_value: &F, x_values: &Vec<F>) -> (Self, LagrangeBasisInnerSteps);
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

#[derive(Debug, Default, Serialize)]
pub struct LagrangeBasisInnerSteps {
    inverse: u128,
}

impl LagrangeBasisInnerSteps {
    fn add_inverse_to_lagrange_basis_inner(&mut self, inverse: u128) {
        self.inverse = inverse;
    }
}

#[derive(Debug, Serialize)]
struct LagrangeBasisAndY {
    x: u128,
    y: u128,
    i: u128,
    lagrange_basis: Vec<u128>,
    steps: LagrangeBasisInnerSteps,
}

#[derive(Debug, Default, Serialize)]
pub struct LagrangeInterpolationSteps {
    y_values: Vec<u128>,
    x_values: Vec<u128>,
    lagrange_basis_and_y_values: Vec<LagrangeBasisAndY>,
    resulting_polynomial: Vec<u128>,
}

impl LagrangeInterpolationSteps {
    fn add_y_values(&mut self, y_values: Vec<u128>) {
        self.y_values = y_values;
    }

    fn add_x_values(&mut self, x_values: Vec<u128>) {
        self.x_values = x_values;
    }

    fn add_to_lagrange_basis_and_y_values(
        &mut self,
        x: u128,
        y: u128,
        i: u128,
        lagrange_basis: Vec<u128>,
        steps: LagrangeBasisInnerSteps,
    ) {
        self.lagrange_basis_and_y_values.push(LagrangeBasisAndY {
            x,
            y,
            lagrange_basis,
            i,
            steps,
        })
    }

    fn add_resulting_polynonial(&mut self, coefficients: Vec<u128>) {
        self.resulting_polynomial = coefficients;
    }
}

impl<
        F: FiniteFieldElement + Clone + Default + Neg<Output = F> + Sub<Output = F> + Add<Output = F>,
    > UnivariatePolynomial<F> for UniPoly<F>
{
    fn coefficients(&self) -> Vec<F> {
        self.coefficients.clone()
    }

    fn evaluate(&self, var: &F) -> F {
        let mut identity = F::zero();
        for (i, x) in self.coefficients.iter().enumerate() {
            let exp = var.pow(i as u128);
            let mul = exp * x.clone();
            identity += mul
        }
        identity
    }

    fn interpolate(y_values: &Vec<F>) -> (Self, LagrangeInterpolationSteps) {
        let mut x_values = vec![];
        for (i, y) in y_values.iter().enumerate() {
            x_values.push(F::new(
                &BigInt::from(i),
                &BigInt::from(y.modulus().unwrap()),
            ));
        }
        Self::interpolate_xy(&x_values, y_values)
    }

    fn interpolate_xy(x_values: &Vec<F>, y_values: &Vec<F>) -> (Self, LagrangeInterpolationSteps) {
        assert_eq!(x_values.len(), y_values.len());
        let mut main_steps = LagrangeInterpolationSteps::default();
        main_steps.add_y_values(
            y_values
                .iter()
                .map(|x| x.element().try_into().unwrap())
                .collect(),
        );
        main_steps.add_x_values(
            x_values
                .iter()
                .map(|x| x.element().try_into().unwrap())
                .collect(),
        );

        let mut resulting_polynomial = Self::zero();
        for (i, (x, y)) in x_values.iter().zip(y_values.iter()).enumerate() {
            let (lagrange_polynomial, steps) = Self::get_lagrange_polynomial(x, &x_values);
            let y_poly = Self::new(vec![y.clone()]);
            let product = &y_poly * &lagrange_polynomial;
            resulting_polynomial = &resulting_polynomial + &product;
            main_steps.add_to_lagrange_basis_and_y_values(
                x.element().try_into().unwrap(),
                y.element().try_into().unwrap(),
                i as u128,
                lagrange_polynomial
                    .coefficients
                    .iter()
                    .map(|x| x.element().try_into().unwrap())
                    .collect(),
                steps,
            );
        }
        main_steps.add_resulting_polynonial(
            resulting_polynomial
                .coefficients()
                .iter()
                .map(|x| x.element().try_into().unwrap())
                .collect(),
        );
        (resulting_polynomial, main_steps)
    }

    fn get_lagrange_polynomial(x_value: &F, x_values: &Vec<F>) -> (Self, LagrangeBasisInnerSteps) {
        /*
         * L_i = \prod_{j=0, j \neq i}^{n} \dfrac{x - x_j}{x_i - x_j}
         *
         *  where:
         *      `i` is x_value
         *      `j` is the index in the loop below
         */
        let mut resulting_polynomial = Self::one();
        let mut running_denominator = F::one();
        let mut lbis = LagrangeBasisInnerSteps::default();
        for x in x_values.iter() {
            if x == x_value {
                continue;
            }
            let numerator = Self::new(vec![-(x.clone()), F::one()]);
            let denominator = x_value.clone() - x.clone();
            running_denominator *= denominator.clone();
            let inverse_of_denominator = Self::new(vec![denominator.inverse().unwrap()]);
            let product = &numerator * &inverse_of_denominator;
            // TODO: implement mul_assign() for &UniPoly
            resulting_polynomial = &resulting_polynomial * &product;
        }
        lbis.add_inverse_to_lagrange_basis_inner(
            running_denominator
                .inverse()
                .unwrap()
                .element()
                .try_into()
                .unwrap(),
        );
        (resulting_polynomial, lbis)
    }
}

impl<F: FiniteFieldElement + Clone> UniPoly<F> {
    pub fn from_latex(poly_string: &str, modulus: &BigInt) -> Self {
        println!("{:?}", poly_string);
        let split_string: Vec<_> = poly_string.split("+").collect();
        let trimmed_split_string: Vec<_> = split_string.iter().map(|x| x.trim()).collect();
        let collections_of_terms_and_coefficients: Vec<_> = trimmed_split_string
            .iter()
            .map(|x| x.split("^").collect::<Vec<_>>())
            .collect();
        let mut degree_and_coefficients = vec![];
        for terms in collections_of_terms_and_coefficients.iter() {
            let exponent: usize;
            let coefficient: usize;
            if terms.len() == 2 {
                let second_characters: Vec<_> = terms[1].chars().collect();
                let exp: usize = second_characters[1].to_string().parse().unwrap();
                exponent = exp;
                if terms[0] == "x" {
                    coefficient = 1;
                } else {
                    let coeff: Vec<_> = terms[0].split("x").collect();
                    coefficient = coeff[0].parse().unwrap();
                }
            } else {
                if terms[0].contains("x") {
                    exponent = 1;
                    let first_characters: Vec<_> = terms[0].chars().collect();
                    if first_characters.len() == 1 {
                        coefficient = 1;
                    } else {
                        coefficient = first_characters[0].to_string().parse().unwrap()
                    }
                } else {
                    exponent = 0;
                    coefficient = terms[0].parse().unwrap();
                }
            }
            degree_and_coefficients.push((exponent, coefficient));
        }
        let degree = degree_and_coefficients[0].0;
        if degree == 0 {
            let coefficients = vec![F::new(&BigInt::from(degree_and_coefficients[0].1), modulus)];
            return UniPoly::new(coefficients);
        } else {
            let mut coefficients = vec![F::zero(); degree + 1];
            for (i, coeff) in coefficients.iter_mut().enumerate() {
                match degree_and_coefficients.iter().find(|x| x.0 == i) {
                    Some(value) => *coeff = F::new(&BigInt::from(value.1), modulus),
                    None => {}
                }
            }
            return UniPoly::new(coefficients);
        }
    }
}

impl<F: FiniteFieldElement + Clone> Mul for &UniPoly<F> {
    type Output = UniPoly<F>;

    fn mul(self, other: Self) -> Self::Output {
        if self.is_zero() || other.is_zero() {
            UniPoly::zero()
        } else {
            let deg_a = self.degree();
            let deg_b = other.degree();
            let max_coefficients = deg_a + deg_b + 1;
            let zero = F::zero();
            let mut product_coefficients = vec![zero; max_coefficients];
            for i in 0..=self.degree() {
                for j in 0..=other.degree() {
                    let index = i + j;
                    let product = self.coefficients[i as usize].clone()
                        * other.coefficients[j as usize].clone();
                    product_coefficients[index as usize] += product;
                }
            }
            let poly = UniPoly::new(product_coefficients);
            poly
        }
    }
}

impl<F: FiniteFieldElement + Clone + Add<Output = F>> Add for &UniPoly<F> {
    type Output = UniPoly<F>;

    fn add(self, rhs: Self) -> Self::Output {
        if self.is_zero() {
            rhs.clone()
        } else if rhs.is_zero() {
            self.clone()
        } else {
            let new_coefficients = add_list(&self.coefficients, &rhs.coefficients);
            UniPoly::new(new_coefficients)
        }
    }
}

fn add_list<T: FiniteFieldElement + Clone + Add<Output = T>>(a: &Vec<T>, b: &Vec<T>) -> Vec<T> {
    let mut res: Vec<T> = vec![];
    if a.len() == b.len() {
        for (x, y) in a.iter().zip(b.iter()) {
            let result = x.clone() + y.clone();
            res.push(result);
        }
    } else if a.len() > b.len() {
        let diff = a.len() - b.len();
        let mut b = b.clone();
        for _ in 0..diff {
            b.push(T::zero());
        }
        for (x, y) in a.iter().zip(b.iter()) {
            let result = x.clone() + y.clone();
            res.push(result);
        }
    } else {
        let diff = b.len() - a.len();
        let mut a = a.clone();
        for _ in 0..diff {
            a.push(T::zero());
        }
        for (x, y) in a.iter().zip(b.iter()) {
            let result = x.clone() + y.clone();
            res.push(result);
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ff::FFE;

    const MODULUS: u128 = 3221225473;

    #[test]
    fn mul() {
        let modulus = BigInt::from(MODULUS);
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
            FFE::new(&BigInt::from(2), &modulus),
            FFE::new(&BigInt::from(-3), &modulus),
            FFE::zero(),
            FFE::new(&BigInt::from(1), &modulus),
        ];
        let poly_1 = UniPoly::<FFE>::new(co_effs_1);
        let co_effs_2 = vec![
            FFE::new(&BigInt::from(5), &modulus),
            FFE::new(&BigInt::from(2), &modulus),
        ];
        let poly_2 = UniPoly::<FFE>::new(co_effs_2);
        let actual = &poly_1 * &poly_2;
        let exp_co_effs = vec![
            FFE::new(&BigInt::from(10), &modulus),
            FFE::new(&BigInt::from(-11), &modulus),
            FFE::new(&BigInt::from(-6), &modulus),
            FFE::new(&BigInt::from(5), &modulus),
            FFE::new(&BigInt::from(2), &modulus),
        ];
        let expected = UniPoly::<FFE>::new(exp_co_effs);
        assert_eq!(actual, expected);
    }

    #[test]
    fn add() {
        let modulus = BigInt::from(MODULUS);
        // Test: (x^2 + x + 5) + (2x^2 + 4x + 2)
        let co_eff_1 = vec![
            FFE::new(&BigInt::from(5), &modulus),
            FFE::new(&BigInt::from(1), &modulus),
            FFE::new(&BigInt::from(1), &modulus),
        ];
        let poly_1 = UniPoly::new(co_eff_1);
        let co_eff_2 = vec![
            FFE::new(&BigInt::from(2), &modulus),
            FFE::new(&BigInt::from(4), &modulus),
            FFE::new(&BigInt::from(2), &modulus),
        ];
        let poly_2 = UniPoly::new(co_eff_2);
        let actual = &poly_1 + &poly_2;
        let exp_co_effs = vec![
            FFE::new(&BigInt::from(7), &modulus),
            FFE::new(&BigInt::from(5), &modulus),
            FFE::new(&BigInt::from(3), &modulus),
        ];
        let expected = UniPoly::<FFE>::new(exp_co_effs);
        assert_eq!(actual, expected);

        // Test: (x^3 - 3x + 2) + (2x + 5)
        let co_eff_3 = vec![
            FFE::new(&BigInt::from(5), &modulus),
            FFE::new(&BigInt::from(2), &modulus),
        ];
        let poly_3 = UniPoly::new(co_eff_3);
        let co_eff_4 = vec![
            FFE::new(&BigInt::from(2), &modulus),
            FFE::new(&BigInt::from(-3), &modulus),
            FFE::new(&BigInt::from(0), &modulus),
            FFE::new(&BigInt::from(1), &modulus),
        ];
        let poly_4 = UniPoly::new(co_eff_4);
        let actual = &poly_3 + &poly_4;
        let exp_co_effs = vec![
            FFE::new(&BigInt::from(7), &modulus),
            FFE::new(&BigInt::from(-1), &modulus),
            FFE::new(&BigInt::from(0), &modulus),
            FFE::new(&BigInt::from(1), &modulus),
        ];
        let expected = UniPoly::<FFE>::new(exp_co_effs);
        assert_eq!(actual, expected);
    }

    #[test]
    fn interpolate() {
        // Interpolating the values: [3, 1, 2, 4]
        let modulus = BigInt::from(17);
        let co_effs = vec![
            FFE::new(&BigInt::from(3), &modulus),
            FFE::new(&BigInt::from(1), &modulus),
            FFE::new(&BigInt::from(2), &modulus),
            FFE::new(&BigInt::from(4), &modulus),
        ];
        let (poly, _steps) = UniPoly::<FFE>::interpolate(&co_effs);
        let exp_co_effs = vec![
            FFE::new(&BigInt::from(3), &modulus),
            FFE::new(&BigInt::from(10), &modulus),
            FFE::new(&BigInt::from(11), &modulus),
            FFE::new(&BigInt::from(11), &modulus),
        ];
        let exp_poly = UniPoly::<FFE>::new(exp_co_effs);
        assert_eq!(exp_poly, poly);
    }

    // #[test]
    // fn test_from_latex() {
    //     let poly_string = "9x^5 + 8x^4 +x^3 + 12x^2 + 3x + 3";
    //     let modulus = BigInt::from(17);
    //     let poly = UniPoly::<FFE>::from_latex(poly_string, &modulus);
    //     println!("{:?}", poly);

    //     let poly_string = "9x^5 + x^3 + 3x + 3";
    //     let modulus = BigInt::from(17);
    //     let poly = UniPoly::<FFE>::from_latex(poly_string, &modulus);
    //     println!("{:?}", poly);

    //     let poly_string = "9x^5 + x^3 + 3x";
    //     let modulus = BigInt::from(17);
    //     let poly = UniPoly::<FFE>::from_latex(poly_string, &modulus);
    //     println!("{:?}", poly);

    //     let poly_string = "13x^4 + 9x^3 + 3x^2 + 8x + 3";
    //     let modulus = BigInt::from(17);
    //     let poly = UniPoly::<FFE>::from_latex(poly_string, &modulus);
    //     println!("{:?}", poly);
    // }
}
