use crate::ff::{FF, FFE};
use crate::uni_poly::{UniPoly, UnivariatePolynomial};

struct FiniteField {
    field: usize,
}

pub fn lagrange_interpolate(
    x_values: &Vec<usize>,
    y_values: &Vec<usize>,
    field: usize,
) -> Vec<usize> {
    let f = field;

    for (x, y) in x_values.iter().zip(y_values.iter()) {}
    todo!()
}
