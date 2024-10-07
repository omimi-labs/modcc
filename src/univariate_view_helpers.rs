use num_bigint::BigInt;

use crate::ff::{FiniteFieldElement, FFE};
use crate::univariate_poly::{LagrangeInterpolationSteps, UniPoly, UnivariatePolynomial};

pub fn evaluate(evaluation_point: isize, poly_string: &str, field: usize) -> usize {
    let modulus = &BigInt::from(field);
    let poly = UniPoly::<FFE>::from_latex(poly_string, modulus);
    let evaluation = poly.evaluate(&FFE::new(&BigInt::from(evaluation_point), modulus));
    return evaluation.element().try_into().unwrap();
}

pub fn lagrange_interpolate(
    x_values: &Vec<i128>,
    y_values: &Vec<i128>,
    field: u128,
) -> (Vec<u128>, LagrangeInterpolationSteps) {
    let mut new_x_values: Vec<FFE> = vec![];
    let mut new_y_values: Vec<FFE> = vec![];
    for (x, y) in x_values.iter().zip(y_values.iter()) {
        let x_ffe = FFE::new(&BigInt::from(*x), &BigInt::from(field));
        new_x_values.push(x_ffe);
        let y_ffe = FFE::new(&BigInt::from(*y), &BigInt::from(field));
        new_y_values.push(y_ffe)
    }
    let (poly, steps) = UniPoly::interpolate_xy(&new_x_values, &new_y_values);
    let coefficients: Vec<u128> = poly
        .coefficients()
        .iter()
        .map(|x| x.element().try_into().unwrap())
        .collect();
    return (coefficients, steps);
}
