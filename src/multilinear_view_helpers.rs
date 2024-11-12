use num_bigint::BigInt;

// use crate::ff::{FiniteFieldElement, FFE};
// use crate::multilinear_poly::{
//     MultilinearLagrangeInterpolationSteps, MultilinearPoly, MultilinearPolynomial, Polynomial,
// };

// pub fn multilinear_interpolate_over_boolean_hypercube(
//     y_values: &Vec<i128>,
//     field: u128,
// ) -> (Vec<(usize, usize)>, MultilinearLagrangeInterpolationSteps) {
//     let mut new_y_values: Vec<FFE> = vec![];
//     for y in y_values.iter() {
//         let y_ffe = FFE::new(&BigInt::from(*y), &BigInt::from(field));
//         new_y_values.push(y_ffe)
//     }
//     let (poly, steps) = MultilinearPoly::interpolate(&new_y_values, field.try_into().unwrap());
//     let coefficients = poly
//         .coefficients()
//         .iter()
//         .map(|(id, x)| (*id, x.element().try_into().unwrap()))
//         .collect::<Vec<(_, _)>>();
//     return (coefficients, steps);
// }
