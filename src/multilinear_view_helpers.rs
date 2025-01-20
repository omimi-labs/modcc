use num_bigint::BigInt;
use num_traits::ToPrimitive;

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
    let mut normalized_eval_len = new_y_values.len().next_power_of_two();
    let num_of_vars: usize;
    if normalized_eval_len == 1 {
        num_of_vars = 1;
        normalized_eval_len = 2;
    } else {
        num_of_vars = normalized_eval_len
            .to_f64()
            .unwrap()
            .log2()
            .to_usize()
            .unwrap();
    }
    let (poly, steps, properties) = create_multilinear_poly(
        &new_y_values,
        &field.try_into().unwrap(),
        num_of_vars,
        normalized_eval_len,
    );
    return (poly.to_latex(num_of_vars), steps, properties);
}
