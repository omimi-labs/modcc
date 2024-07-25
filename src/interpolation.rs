use crate::ff::{FiniteFieldElement, FFE};
use crate::uni_poly::{UniPoly, UnivariatePolynomial};

pub fn lagrange_interpolate(
    x_values: &Vec<isize>,
    y_values: &Vec<isize>,
    field: usize,
) -> Vec<usize> {
    let mut new_x_values: Vec<FFE> = vec![];
    let mut new_y_values: Vec<FFE> = vec![];
    for (x, y) in x_values.iter().zip(y_values.iter()) {
        let x_ffe = FFE::new(*x, field);
        new_x_values.push(x_ffe);
        let y_ffe = FFE::new(*y, field);
        new_y_values.push(y_ffe)
    }
    let poly = UniPoly::interpolate_xy(&new_x_values, &new_y_values);
    let mut coefficients = vec![];
    for x in poly.coefficients().iter() {
        coefficients.push(x.element())
    }
    return coefficients;
}
