use std::cmp::max;
use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::ops::{Add, Mul, Neg};
use std::{u128, usize};

use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};
use serde::Serialize;

use crate::ff::FiniteFieldElement;

pub trait Polynomial<F>: Sized {
    fn is_zero(&self) -> bool;

    fn is_one(&self) -> bool;

    fn new(num_of_vars: usize, coefficients: BTreeMap<usize, F>) -> Self;

    fn zero() -> Self;

    fn one() -> Self;

    fn coefficients(&self) -> &BTreeMap<usize, F>;
}

#[derive(Debug, Default, Clone, Serialize)]
pub struct MultilinearLagrangeInterpolationSteps {
    y_values: Vec<u128>,
    binary_strings: Vec<String>,
    collection_of_ids_and_coefficients: Vec<Vec<(usize, u128)>>,
}

impl MultilinearLagrangeInterpolationSteps {
    fn add_y_values(&mut self, y_values: Vec<u128>) {
        self.y_values = y_values;
    }

    fn add_binary_strings(&mut self, binary_string: String) {
        self.binary_strings.push(binary_string);
    }

    fn add_coefficients(&mut self, id_and_coefficients: Vec<(usize, u128)>) {
        self.collection_of_ids_and_coefficients
            .push(id_and_coefficients);
    }
}

pub trait MultilinearPolynomial<F>: Polynomial<F> {
    fn scalar_mul(&mut self, scalar: F);
    fn interpolate(
        y_values: &Vec<F>,
        field: BigInt,
    ) -> (Self, MultilinearLagrangeInterpolationSteps);
}

#[derive(Debug, Clone, PartialEq)]
pub struct MultilinearPoly<F> {
    num_of_vars: usize,
    coefficients: BTreeMap<usize, F>,
}

impl<F: FiniteFieldElement> Polynomial<F> for MultilinearPoly<F> {
    fn new(num_of_vars: usize, coefficients: BTreeMap<usize, F>) -> Self {
        Self {
            num_of_vars,
            coefficients,
        }
    }

    fn is_zero(&self) -> bool {
        return self.num_of_vars == 0 && self.coefficients.is_empty();
    }

    fn is_one(&self) -> bool {
        return self.num_of_vars == 0 && !self.coefficients.is_empty();
    }

    fn zero() -> Self {
        let coefficients = BTreeMap::new();
        Self::new(0, coefficients)
    }

    fn one() -> Self {
        let mut coefficients = BTreeMap::new();
        coefficients.insert(0, F::one());
        Self::new(0, coefficients)
    }

    fn coefficients(&self) -> &BTreeMap<usize, F> {
        &self.coefficients
    }
}

impl<F: FiniteFieldElement + Clone + Neg<Output = F> + Add<Output = F>> MultilinearPolynomial<F>
    for MultilinearPoly<F>
{
    fn scalar_mul(&mut self, scalar: F) {
        for value in self.coefficients.values_mut() {
            *value *= scalar.clone();
        }
    }

    fn interpolate(
        y_values: &Vec<F>,
        field: BigInt,
    ) -> (Self, MultilinearLagrangeInterpolationSteps) {
        let next_power_of_two = y_values.len().next_power_of_two();
        let bit_size = next_power_of_two
            .to_f64()
            .unwrap()
            .log2()
            .to_usize()
            .unwrap();
        let mut evaluations = y_values.clone();
        evaluations.resize(next_power_of_two, F::zero());
        let mut multilinear_interpolation_steps = MultilinearLagrangeInterpolationSteps::default();
        multilinear_interpolation_steps.add_y_values(
            evaluations
                .iter()
                .map(|x| x.element().try_into().unwrap())
                .collect(),
        );
        let bit_iterator = LagrangeBasisBooleanHyperCube::<F>::new(bit_size, &field);
        let mut result = Self::zero();
        for ((vec_of_poly, binary_string), y_value) in bit_iterator.zip(evaluations.iter()) {
            multilinear_interpolation_steps.add_binary_strings(binary_string);
            let mut inter_poly = Self::one();
            for poly in vec_of_poly.iter() {
                let res = &inter_poly * poly;
                inter_poly = res;
            }
            multilinear_interpolation_steps.add_coefficients(
                inter_poly
                    .coefficients
                    .iter()
                    .map(|(key, value)| (key.clone(), value.element().try_into().unwrap()))
                    .collect::<Vec<(_, _)>>(),
            );
            inter_poly.scalar_mul(y_value.clone());
            result = &result + &inter_poly;
        }

        (result, multilinear_interpolation_steps)
    }
}

impl<F: FiniteFieldElement + Clone> Mul for &MultilinearPoly<F> {
    type Output = MultilinearPoly<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.is_zero() || rhs.is_zero() {
            return MultilinearPoly::zero();
        }
        if self.is_one() {
            return rhs.clone();
        }
        if rhs.is_one() {
            return self.clone();
        }
        let mut new_coefficients = BTreeMap::new();
        for (i, self_co_eff) in self.coefficients.iter() {
            for (j, rhs_co_eff) in rhs.coefficients.iter() {
                let i_bits = binary_string(*i, self.num_of_vars);
                let j_bits = binary_string(*j, rhs.num_of_vars);
                let combined_bits = i_bits + &j_bits;
                let index = usize::from_str_radix(&combined_bits, 2).unwrap();
                let co_eff = self_co_eff.clone() * rhs_co_eff.clone();
                new_coefficients.insert(index, co_eff);
            }
        }
        let new_multi_poly =
            MultilinearPoly::new(self.num_of_vars + rhs.num_of_vars, new_coefficients);
        new_multi_poly
    }
}

impl<F: FiniteFieldElement + Clone + Add<Output = F>> Add for &MultilinearPoly<F> {
    type Output = MultilinearPoly<F>;

    fn add(self, rhs: Self) -> Self::Output {
        if self.is_zero() {
            return rhs.clone();
        } else if rhs.is_zero() {
            return self.clone();
        } else {
            let new_coefficients = add_coefficients(
                &self.coefficients,
                self.num_of_vars,
                &rhs.coefficients,
                rhs.num_of_vars,
            );
            let new_multi_poly =
                MultilinearPoly::new(max(self.num_of_vars, rhs.num_of_vars), new_coefficients);
            new_multi_poly
        }
    }
}

fn add_coefficients<T: FiniteFieldElement + Clone + Add<Output = T>>(
    a: &BTreeMap<usize, T>,
    n: usize,
    b: &BTreeMap<usize, T>,
    m: usize,
) -> BTreeMap<usize, T> {
    let mut new_coefficients = BTreeMap::new();
    if n == m {
        for ((i, a_coeff), (_, b_coeff)) in a.iter().zip(b.iter()) {
            new_coefficients.insert(*i, a_coeff.clone() + b_coeff.clone());
        }
    } else if n > m {
        let zero = T::zero();
        for (i, a_coeff) in a.iter() {
            let b_coeff = b.get(i).unwrap_or(&zero);
            new_coefficients.insert(*i, a_coeff.clone() + b_coeff.clone());
        }
    } else {
        let zero = T::zero();
        for (i, b_coeff) in b.iter() {
            let a_coeff = b.get(i).unwrap_or(&zero);
            new_coefficients.insert(*i, a_coeff.clone() + b_coeff.clone());
        }
    }
    new_coefficients
}

/// Convert a number to a binary string of a given size
pub fn binary_string(value: usize, bit_count: usize) -> String {
    let binary = format!("{:b}", value);
    "0".repeat(bit_count.saturating_sub(binary.len())) + &binary
}

pub struct LagrangeBasisBooleanHyperCube<F> {
    bit_size: usize,
    total_points: usize,
    current_point: usize,
    modulus: BigInt,
    _marker: PhantomData<F>,
}

impl<F: FiniteFieldElement> LagrangeBasisBooleanHyperCube<F> {
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

impl<F: FiniteFieldElement + Clone + Neg<Output = F>> Iterator
    for LagrangeBasisBooleanHyperCube<F>
{
    type Item = (Vec<MultilinearPoly<F>>, String);

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
            .map(|a| {
                if a == '1' {
                    let mut btree_map = BTreeMap::new();
                    btree_map.insert(0, F::new(&zero, &self.modulus));
                    btree_map.insert(1, F::new(&one, &self.modulus));
                    MultilinearPoly::new(1, btree_map)
                } else {
                    let mut btree_map = BTreeMap::new();
                    let minus_one = -F::new(&one, &self.modulus);
                    btree_map.insert(0, F::new(&one, &self.modulus));
                    btree_map.insert(1, minus_one);
                    MultilinearPoly::new(1, btree_map)
                }
            })
            .collect::<Vec<MultilinearPoly<F>>>();

        self.current_point += 1;

        Some((point, index_as_binary))
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use num_bigint::BigInt;
    use num_traits::{One, Zero};

    use super::{MultilinearPoly, MultilinearPolynomial, Polynomial};
    use crate::ff::{FiniteFieldElement, FFE};
    use crate::multilinear_poly::LagrangeBasisBooleanHyperCube;

    #[test]
    fn test_lagrange_basis_boolean_hypercube_iteration() {
        let modulus = BigInt::from(17);

        let mut two_bit_iterator = LagrangeBasisBooleanHyperCube::<FFE>::new(2, &modulus);
        assert_eq!(two_bit_iterator.total_points, 4);

        let zero = BigInt::zero();
        let one = BigInt::one();
        let mut btree_map = BTreeMap::new();
        let minus_one = -FFE::new(&one, &modulus);
        btree_map.insert(0, FFE::new(&one, &modulus));
        btree_map.insert(1, minus_one);
        let check_zero_poly = MultilinearPoly::new(1, btree_map);

        let mut btree_map = BTreeMap::new();
        btree_map.insert(0, FFE::new(&zero, &modulus));
        btree_map.insert(1, FFE::new(&one, &modulus));
        let check_one_poly = MultilinearPoly::new(1, btree_map);
        assert_eq!(
            two_bit_iterator.next(),
            Some((
                vec![check_zero_poly.clone(), check_zero_poly.clone()],
                String::from("00")
            ))
        );
        assert_eq!(
            two_bit_iterator.next(),
            Some((
                vec![check_zero_poly.clone(), check_one_poly.clone()],
                String::from("01")
            ))
        );
        assert_eq!(
            two_bit_iterator.next(),
            Some((
                vec![check_one_poly.clone(), check_zero_poly.clone()],
                String::from("10")
            ))
        );
        assert_eq!(
            two_bit_iterator.next(),
            Some((
                vec![check_one_poly.clone(), check_one_poly.clone()],
                String::from("11")
            ))
        );
        assert_eq!(two_bit_iterator.next(), None);

        let mut three_bit_iterator = LagrangeBasisBooleanHyperCube::<FFE>::new(3, &modulus);
        assert_eq!(three_bit_iterator.total_points, 8);
        assert_eq!(
            three_bit_iterator.next(),
            Some((
                vec![
                    check_zero_poly.clone(),
                    check_zero_poly.clone(),
                    check_zero_poly.clone()
                ],
                String::from("000")
            ))
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some((
                vec![
                    check_zero_poly.clone(),
                    check_zero_poly.clone(),
                    check_one_poly.clone()
                ],
                String::from("001")
            ))
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some((
                vec![
                    check_zero_poly.clone(),
                    check_one_poly.clone(),
                    check_zero_poly.clone()
                ],
                String::from("010")
            ))
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some((
                vec![
                    check_zero_poly.clone(),
                    check_one_poly.clone(),
                    check_one_poly.clone()
                ],
                String::from("011")
            ))
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some((
                vec![
                    check_one_poly.clone(),
                    check_zero_poly.clone(),
                    check_zero_poly.clone()
                ],
                String::from("100")
            ))
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some((
                vec![
                    check_one_poly.clone(),
                    check_zero_poly.clone(),
                    check_one_poly.clone()
                ],
                String::from("101")
            ))
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some((
                vec![
                    check_one_poly.clone(),
                    check_one_poly.clone(),
                    check_zero_poly
                ],
                String::from("110")
            ))
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some((
                vec![
                    check_one_poly.clone(),
                    check_one_poly.clone(),
                    check_one_poly
                ],
                String::from("111")
            ))
        );
        assert_eq!(three_bit_iterator.next(), None);
    }

    #[test]
    fn test_interpolate() {
        let modulus = BigInt::from(17);
        // x + y + z + 7
        let evaluations = vec![
            FFE::new(&BigInt::from(7), &modulus),
            FFE::new(&BigInt::from(8), &modulus),
            FFE::new(&BigInt::from(8), &modulus),
            FFE::new(&BigInt::from(9), &modulus),
            FFE::new(&BigInt::from(8), &modulus),
            FFE::new(&BigInt::from(9), &modulus),
            FFE::new(&BigInt::from(9), &modulus),
            FFE::new(&BigInt::from(10), &modulus),
        ];
        let (poly, _) = MultilinearPoly::interpolate(&evaluations, modulus.clone());

        let num_of_vars = 3;
        let mut coefficients = BTreeMap::new();
        coefficients.insert(0, FFE::new(&BigInt::from(7), &modulus));
        coefficients.insert(1, FFE::new(&BigInt::from(1), &modulus));
        coefficients.insert(2, FFE::new(&BigInt::from(1), &modulus));
        coefficients.insert(4, FFE::new(&BigInt::from(1), &modulus));
        coefficients.insert(3, FFE::new(&BigInt::from(0), &modulus));
        coefficients.insert(5, FFE::new(&BigInt::from(0), &modulus));
        coefficients.insert(6, FFE::new(&BigInt::from(0), &modulus));
        coefficients.insert(7, FFE::new(&BigInt::from(0), &modulus));

        let expected_poly = MultilinearPoly::new(num_of_vars, coefficients);
        assert_eq!(poly, expected_poly);
    }
}
