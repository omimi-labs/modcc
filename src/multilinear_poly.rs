use std::cmp::max;
use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::ops::{Add, Mul, Neg};

use num_traits::ToPrimitive;

use crate::ff::FiniteFieldElement;

pub trait Polynomial<F>: Sized {
    fn is_zero(&self) -> bool;

    fn is_one(&self) -> bool;

    fn new(num_of_vars: usize, coefficients: BTreeMap<usize, F>) -> Self;

    fn zero() -> Self;

    fn one() -> Self;
}

pub trait MultilinearPolynomial<F>: Polynomial<F> {
    fn interpolate(y_values: &Vec<F>) -> Self;
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
}

impl<F: FiniteFieldElement + Clone + Neg<Output = F> + Add<Output = F>> MultilinearPolynomial<F>
    for MultilinearPoly<F>
{
    fn interpolate(y_values: &Vec<F>) -> Self {
        let next_power_of_two = y_values.len().next_power_of_two();
        let bit_size = next_power_of_two
            .to_f64()
            .unwrap()
            .log2()
            .to_usize()
            .unwrap();
        let mut bit_iterator = LagrangeBasisBooleanHyperCube::<F>::new(bit_size);
        let result = Self::zero();
        for vec_of_poly in bit_iterator {
            let mut inter_poly = Self::one();
            for poly in vec_of_poly.iter() {
                let res = &inter_poly + poly;
                inter_poly = res;
            }
            // result = result + inter_poly;
        }
        todo!()
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
    _marker: PhantomData<F>,
}

impl<F: FiniteFieldElement> LagrangeBasisBooleanHyperCube<F> {
    pub fn new(bit_size: usize) -> Self {
        Self {
            bit_size,
            total_points: 2_usize.pow(bit_size as u32),
            current_point: 0,
            _marker: PhantomData,
        }
    }
}

impl<F: FiniteFieldElement + Clone + Neg<Output = F>> Iterator
    for LagrangeBasisBooleanHyperCube<F>
{
    type Item = Vec<MultilinearPoly<F>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_point == self.total_points || self.bit_size == 0 {
            return None;
        }

        // convert the current index to binary value of the given length
        let index_as_binary = binary_string(self.current_point, self.bit_size);
        let point = index_as_binary
            .chars()
            .map(|a| {
                if a == '1' {
                    let mut btree_map = BTreeMap::new();
                    btree_map.insert(0, F::zero());
                    btree_map.insert(1, F::one());
                    MultilinearPoly::new(1, btree_map)
                } else {
                    let mut btree_map = BTreeMap::new();
                    let minus_one = -F::one();
                    btree_map.insert(0, F::one());
                    btree_map.insert(1, minus_one);
                    MultilinearPoly::new(1, btree_map)
                }
            })
            .collect::<Vec<MultilinearPoly<F>>>();

        self.current_point += 1;

        Some(point)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::{MultilinearPoly, Polynomial};
    use crate::ff::{FiniteFieldElement, FFE};
    use crate::multilinear_poly::LagrangeBasisBooleanHyperCube;

    #[test]
    fn test_lagrange_basis_boolean_hypercube_iteration() {
        let mut two_bit_iterator = LagrangeBasisBooleanHyperCube::<FFE>::new(2);
        assert_eq!(two_bit_iterator.total_points, 4);

        let mut btree_map = BTreeMap::new();
        let minus_one = -FFE::one();
        btree_map.insert(0, FFE::one());
        btree_map.insert(1, minus_one);
        let check_zero_poly = MultilinearPoly::new(1, btree_map);

        let mut btree_map = BTreeMap::new();
        btree_map.insert(0, FFE::zero());
        btree_map.insert(1, FFE::one());
        let check_one_poly = MultilinearPoly::new(1, btree_map);
        assert_eq!(
            two_bit_iterator.next(),
            Some(vec![check_zero_poly.clone(), check_zero_poly.clone()])
        );
        assert_eq!(
            two_bit_iterator.next(),
            Some(vec![check_zero_poly.clone(), check_one_poly.clone()])
        );
        assert_eq!(
            two_bit_iterator.next(),
            Some(vec![check_one_poly.clone(), check_zero_poly.clone()])
        );
        assert_eq!(
            two_bit_iterator.next(),
            Some(vec![check_one_poly.clone(), check_one_poly.clone()])
        );
        assert_eq!(two_bit_iterator.next(), None);

        let mut three_bit_iterator = LagrangeBasisBooleanHyperCube::<FFE>::new(3);
        assert_eq!(three_bit_iterator.total_points, 8);
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![
                check_zero_poly.clone(),
                check_zero_poly.clone(),
                check_zero_poly.clone()
            ])
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![
                check_zero_poly.clone(),
                check_zero_poly.clone(),
                check_one_poly.clone()
            ])
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![
                check_zero_poly.clone(),
                check_one_poly.clone(),
                check_zero_poly.clone()
            ])
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![
                check_zero_poly.clone(),
                check_one_poly.clone(),
                check_one_poly.clone()
            ])
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![
                check_one_poly.clone(),
                check_zero_poly.clone(),
                check_zero_poly.clone()
            ])
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![
                check_one_poly.clone(),
                check_zero_poly.clone(),
                check_one_poly.clone()
            ])
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![
                check_one_poly.clone(),
                check_one_poly.clone(),
                check_zero_poly
            ])
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![
                check_one_poly.clone(),
                check_one_poly.clone(),
                check_one_poly
            ])
        );
        assert_eq!(three_bit_iterator.next(), None);
    }

    #[test]
    fn test_mul() {
        let mut check_one_coeffs = BTreeMap::new();
        check_one_coeffs.insert(0, FFE::zero());
        check_one_coeffs.insert(1, FFE::one());
        let check_one = MultilinearPoly::new(1, check_one_coeffs);

        let mut check_zero_coeffs = BTreeMap::new();
        let minus_one = -FFE::one();
        check_zero_coeffs.insert(0, FFE::one());
        check_zero_coeffs.insert(1, minus_one);
        let check_zero = MultilinearPoly::new(1, check_zero_coeffs);

        let multi_poly_1 = &check_zero * &check_zero;

        let multi_poly_2 = &check_zero * &check_one;

        let multi_poly_3 = &check_one * &check_zero;

        let multi_poly_4 = &check_one * &check_one;

        let multi_poly_5 = &check_zero * &multi_poly_1;
    }
}
