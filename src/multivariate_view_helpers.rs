use num_bigint::BigInt;

use crate::ff::{FiniteFieldElement, FFE};
use crate::multivariate_poly::{MultivariatePoly, MultivariatePolynomial, Properties, Steps};

pub fn multivariate_interpolate_over_finite_field(
    num_of_vars: usize,
    evaluation_points: &Vec<Vec<i128>>,
    y_values: &Vec<i128>,
    field: usize,
) -> Result<(String, Steps, Properties), String> {
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
    let (poly, steps, properties) =
        MultivariatePoly::interpolate(&group_of_evaluation_points, &new_y_values);
    Ok((poly.to_latex(num_of_vars), steps, properties))
}

pub fn multivariate_full_evaluation(
    evaluation_points: &Vec<(usize, usize)>,
    poly_string: &String,
    field: usize,
) -> String {
    let poly = MultivariatePoly::<FFE>::from_latex(poly_string.clone(), &BigInt::from(field));
    let eval_points = evaluation_points
        .iter()
        .map(|(i, val)| (*i, FFE::new(&BigInt::from(*val), &BigInt::from(field))))
        .collect::<Vec<(usize, FFE)>>();
    let evaluation = poly.full_evaluation(&eval_points);
    format!("${}$", evaluation.element())
}

pub fn multivariate_partial_evaluation(
    evaluation_points: &Vec<(usize, usize)>,
    poly_string: &String,
    field: usize,
) -> String {
    let poly = MultivariatePoly::<FFE>::from_latex(poly_string.clone(), &BigInt::from(field));
    let eval_points = evaluation_points
        .iter()
        .map(|(i, val)| (*i, FFE::new(&BigInt::from(*val), &BigInt::from(field))))
        .collect::<Vec<(usize, FFE)>>();
    let partial_evaluation = poly.partial_evaluation(&eval_points);
    format!("${}$", partial_evaluation.create_latex_rhs())
}
