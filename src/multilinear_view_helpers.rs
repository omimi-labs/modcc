use num_bigint::BigInt;

use crate::ff::{FiniteFieldElement, FFE};
use crate::multilinear_poly::create_multilinear_poly;
use crate::multivariate_poly::{MultivariatePolynomial, Properties, Steps};

pub fn multilinear_interpolate_over_boolean_hypercube(
    y_values: &Vec<i128>,
    field: u128,
) -> (String, Steps, Properties) {
    let mut new_y_values = vec![];
    for y in y_values.iter() {
        let y_ffe = FFE::new(&BigInt::from(*y), &BigInt::from(field));
        new_y_values.push(y_ffe)
    }
    let (poly, steps, properties) =
        create_multilinear_poly(&new_y_values, &field.try_into().unwrap());
    return (poly.to_latex(), steps, properties);
}
