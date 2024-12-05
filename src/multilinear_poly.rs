use std::marker::PhantomData;
use std::ops::{Add, Neg, Sub};
use std::usize;

use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};

use crate::ff::FiniteFieldElement;
use crate::multivariate_poly::{MultivariatePoly, MultivariatePolynomial, Properties, Steps};

/// Convert a number to a binary string of a given size
pub fn binary_string(value: usize, bit_count: usize) -> String {
    let binary = format!("{:b}", value);
    "0".repeat(bit_count.saturating_sub(binary.len())) + &binary
}

pub struct BooleanHyperCube<F> {
    bit_size: usize,
    total_points: usize,
    current_point: usize,
    modulus: BigInt,
    _marker: PhantomData<F>,
}

impl<F: FiniteFieldElement> BooleanHyperCube<F> {
    pub fn new(bit_size: usize, field: &BigInt) -> Self {
        Self {
            bit_size,
            total_points: 2_usize.pow(bit_size as u32),
            current_point: 0,
            modulus: field.clone(),
            _marker: PhantomData,
        }
    }
}

impl<F: FiniteFieldElement + Clone + Neg<Output = F>> Iterator for BooleanHyperCube<F> {
    type Item = (Vec<F>, String);

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_point == self.total_points || self.bit_size == 0 {
            return None;
        }

        // convert the current index to binary value of the given length
        let index_as_binary = binary_string(self.current_point, self.bit_size);
        let zero = BigInt::zero();
        let one = BigInt::one();

        let point = index_as_binary
            .chars()
            .map(|x| {
                if x == '1' {
                    F::new(&one, &self.modulus)
                } else {
                    F::new(&zero, &self.modulus)
                }
            })
            .collect::<Vec<F>>();
        self.current_point += 1;
        Some((point, index_as_binary))
    }
}

pub fn create_multilinear_poly<
    F: FiniteFieldElement + Neg<Output = F> + Sub<Output = F> + Add<Output = F> + Clone,
>(
    y_values: &Vec<F>,
    field: &BigInt,
) -> (MultivariatePoly<F>, Steps, Properties) {
    let normalized_eval_len = y_values.len().next_power_of_two();
    let num_of_vars = normalized_eval_len
        .to_f64()
        .unwrap()
        .log2()
        .to_usize()
        .unwrap();
    let mut evaluations = y_values.clone();
    evaluations.resize(normalized_eval_len, F::zero());
    let bit_iterator = BooleanHyperCube::<F>::new(num_of_vars, field);
    let mut group_of_evaluation_points = vec![];
    for (evaluation_points, _) in bit_iterator {
        group_of_evaluation_points.push(evaluation_points);
    }
    MultivariatePoly::interpolate(&group_of_evaluation_points, &evaluations)
}
