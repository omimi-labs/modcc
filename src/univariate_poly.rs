use std::{
    fmt::Debug,
    ops::{Add, Mul, Neg, Sub},
    u128,
};

use crate::ff::FiniteFieldElement;
use num_bigint::BigInt;
use num_traits::{One, Zero};
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
#[derive(PartialEq, Clone)]
pub struct UniPoly<F> {
    // Co-efficient represented from lower degree to higher
    // For example: 2x^2 + x + 1 is represented as [1, 1, 2]
    coefficients: Vec<F>,
}

impl<F: FiniteFieldElement + Clone + Add<Output = F>> Debug for UniPoly<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_latex())?;
        Ok(())
    }
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
        if x_values.len() == 1 {
            // TODO: Add steps for when this case
            return (Self::new(y_values.to_vec()), main_steps);
        }
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
        let poly_string_rhs: Vec<_> = poly_string.trim_matches('$').split("=").collect();
        let split_string: Vec<_> = poly_string_rhs[1].split("+").collect();
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
                let exp: usize = terms[1].parse().unwrap();
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
                    if terms[0].trim_end_matches("x") == "" {
                        coefficient = 1;
                    } else {
                        coefficient = terms[0].trim_end_matches("x").parse().unwrap()
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

    pub fn to_latex(&self) -> String {
        let mut latex_string = String::from("f(x) = ");
        for (i, value) in self.coefficients.iter().rev().enumerate() {
            if value.is_zero() {
                if self.coefficients.len() == 1 {
                    latex_string += "0"
                }
            } else if i == self.coefficients.len() - 1 {
                latex_string += &format!("{:?}", value);
            } else {
                if value.is_one() {
                    if i == self.coefficients.len() - 2 {
                        latex_string += &format!("x + ")
                    } else {
                        latex_string += &format!("x^{} + ", self.coefficients.len() - 1 - i)
                    }
                } else {
                    if i == self.coefficients.len() - 2 {
                        latex_string += &format!("{:?}x + ", value)
                    } else {
                        latex_string +=
                            &format!("{:?}x^{} + ", value, self.coefficients.len() - 1 - i)
                    }
                }
            }
        }
        format!("${}$", latex_string.trim().trim_end_matches("+").trim())
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

    use rand::{thread_rng, Rng};

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

    fn generate_x_points(num_of_points: usize, modulus: &BigInt) -> Vec<FFE> {
        let mut x_points = vec![];
        let mut rng = thread_rng();
        while x_points.len() != num_of_points {
            let x_point: usize = rng.gen();
            let element = BigInt::from(x_point);
            let var = FFE::new(&element, modulus);
            if x_points.contains(&var) {
                continue;
            }
            x_points.push(var);
        }
        x_points
    }

    fn generate_y_points(num_of_points: usize, modulus: &BigInt) -> Vec<FFE> {
        let mut y_points = vec![];
        let mut rng = thread_rng();
        while y_points.len() != num_of_points {
            let y_point: usize = rng.gen();
            let element = BigInt::from(y_point);
            let var = FFE::new(&element, modulus);
            y_points.push(var);
        }
        y_points
    }

    #[test]
    fn interpolate() {
        let modulus = BigInt::from(17);
        let rounds = 10;
        for num_of_points in 1..100 {
            for _ in 0..rounds {
                let x_values = generate_x_points(num_of_points.try_into().unwrap(), &modulus);
                let y_values = generate_y_points(num_of_points.try_into().unwrap(), &modulus);
                let (poly, _) = UniPoly::interpolate_xy(&x_values, &y_values);
                for (x_value, y_value) in x_values.iter().zip(y_values.iter()) {
                    let evaluation = poly.evaluate(&x_value);
                    assert_eq!(evaluation, *y_value)
                }
            }
        }
    }

    #[test]
    fn test_latex() {
        let modulus = BigInt::from(17);
        let rounds = 10;
        for num_of_points in 1..100 {
            for _ in 0..rounds {
                let x_values = generate_x_points(num_of_points.try_into().unwrap(), &modulus);
                let y_values = generate_y_points(num_of_points.try_into().unwrap(), &modulus);
                let (poly, _) = UniPoly::interpolate_xy(&x_values, &y_values);
                let poly_latex_string = poly.to_latex();
                let gen_poly = UniPoly::<FFE>::from_latex(&poly_latex_string, &modulus);
                for (x_value, y_value) in x_values.iter().zip(y_values.iter()) {
                    let poly_evaluation = poly.evaluate(&x_value);
                    let gen_poly_evaluation = gen_poly.evaluate(&x_value);
                    assert_eq!(poly_evaluation, *y_value);
                    assert_eq!(gen_poly_evaluation, *y_value);
                }
            }
        }
    }
}
