use num_bigint::BigInt;

use crate::ff::{FiniteFieldElement, FFE};
use crate::multivariate_poly::{MultivariatePoly, MultivariatePolynomial, Polynomial};

pub fn multivariate_interpolate_over_finite_field(
    num_of_vars: usize,
    evaluation_points: &Vec<Vec<i128>>,
    y_values: &Vec<i128>,
    field: usize,
) -> Result<Vec<(usize, Vec<(usize, usize)>)>, String> {
    for group in evaluation_points {
        if group.len() != num_of_vars {
            return Err(String::from("Invalid collection of vars"));
        }
    }
    if evaluation_points.len() != y_values.len() {
        return Err(String::from("Invalid length"));
    }
    let modulus = BigInt::from(field);
    let new_y_values = y_values
        .iter()
        .map(|x| FFE::new(&BigInt::from(x.clone()), &modulus))
        .collect::<Vec<FFE>>();
    let group_of_evaluation_points: Vec<Vec<FFE>> = evaluation_points
        .iter()
        .map(|x| {
            {
                x.iter()
                    .map(|y| FFE::new(&BigInt::from(y.clone()), &modulus))
            }
            .collect::<Vec<_>>()
        })
        .collect();
    let poly = MultivariatePoly::interpolate(&group_of_evaluation_points, &new_y_values);
    println!("{:?}", poly);
    let terms = poly
        .terms()
        .iter()
        .map(|term| {
            (
                term.coefficient().element().try_into().unwrap(),
                term.vars()
                    .iter()
                    .map(|x| (x.var_index(), x.power()))
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<Vec<(usize, Vec<_>)>>();
    return Ok(terms);
}
