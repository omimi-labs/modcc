use std::marker::PhantomData;
use std::ops::Neg;

use crate::ff::FiniteFieldElement;
use crate::univariate_poly::{Polynomial as Poly, UniPoly};

/// Structure for point iteration over boolean hypercube
/// e.g. BooleanHyperCube 2 variables
/// Some(00), Some(01), Some(10), Some(11), None
pub struct BooleanHyperCube<F> {
    bit_size: usize,
    total_points: usize,
    current_point: usize,
    _marker: PhantomData<F>,
}

impl<F: FiniteFieldElement> BooleanHyperCube<F> {
    pub fn new(bit_size: usize) -> Self {
        Self {
            bit_size,
            total_points: 2_usize.pow(bit_size as u32),
            current_point: 0,
            _marker: PhantomData,
        }
    }
}

/// Convert a number to a binary string of a given size
pub fn binary_string(index: usize, bit_count: usize) -> String {
    let binary = format!("{:b}", index);
    "0".repeat(bit_count.saturating_sub(binary.len())) + &binary
}

impl<F: FiniteFieldElement> Iterator for BooleanHyperCube<F> {
    type Item = Vec<F>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_point == self.total_points || self.bit_size == 0 {
            return None;
        }

        // convert the current index to binary value of the given length
        let index_as_binary = binary_string(self.current_point, self.bit_size);
        let point = index_as_binary
            .chars()
            .map(|a| if a == '1' { F::one() } else { F::zero() })
            .collect::<Vec<F>>();

        self.current_point += 1;

        Some(point)
    }
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
    type Item = Vec<UniPoly<F>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_point == self.total_points || self.bit_size == 0 {
            return None;
        }

        // convert the current index to binary value of the given length
        let index_as_binary = binary_string(self.current_point, self.bit_size);
        let minus_one = -F::one();
        let check_zero_poly = UniPoly::<F>::new(vec![F::one(), minus_one]);
        let check_one_poly = UniPoly::<F>::new(vec![F::zero(), F::one()]);
        let point = index_as_binary
            .chars()
            .map(|a| {
                if a == '1' {
                    check_one_poly.clone()
                } else {
                    check_zero_poly.clone()
                }
            })
            .collect::<Vec<UniPoly<F>>>();

        self.current_point += 1;

        Some(point)
    }
}

pub trait Polynomial<F>: Sized {
    fn is_zero(&self) -> bool;

    fn degree(&self) -> u128;

    fn new(coefficients: Vec<F>) -> Self;

    fn x() -> Self;

    fn zero() -> Self;

    fn one() -> Self;
}

pub trait MultilinearPolynomial<F: Default>: Polynomial<F> {
    fn evaluate(&self, x: F) -> F;

    fn interpolate(y_values: &Vec<F>) -> Self;
}

struct LagrangeBasis<F> {
    values: Vec<UniPoly<F>>,
}

struct MultilinearPoly<F> {
    num_of_vars: u128,
    _dummy: PhantomData<F>,
}

impl<F: FiniteFieldElement + Neg<Output = F>> MultilinearPoly<F> {
    fn _interpolate(_num_of_vars: u128, _y_values: &Vec<F>) -> Self {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use crate::ff::{FiniteFieldElement, FFE};
    use crate::multilinear_poly::{BooleanHyperCube, LagrangeBasisBooleanHyperCube};
    use crate::univariate_poly::{Polynomial, UniPoly};

    #[test]
    fn test_boolean_hypercube_iteration() {
        let mut two_bit_iterator = BooleanHyperCube::<FFE>::new(2);
        assert_eq!(two_bit_iterator.total_points, 4);
        assert_eq!(two_bit_iterator.next(), Some(vec![FFE::zero(); 2]));
        assert_eq!(two_bit_iterator.next(), Some(vec![FFE::zero(), FFE::one()]));
        assert_eq!(two_bit_iterator.next(), Some(vec![FFE::one(), FFE::zero()]));
        assert_eq!(two_bit_iterator.next(), Some(vec![FFE::one(); 2]));
        assert_eq!(two_bit_iterator.next(), None);

        let mut three_bit_iterator = BooleanHyperCube::<FFE>::new(3);
        assert_eq!(three_bit_iterator.total_points, 8);
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![FFE::zero(), FFE::zero(), FFE::zero()])
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![FFE::zero(), FFE::zero(), FFE::one()])
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![FFE::zero(), FFE::one(), FFE::zero()])
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![FFE::zero(), FFE::one(), FFE::one()])
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![FFE::one(), FFE::zero(), FFE::zero()])
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![FFE::one(), FFE::zero(), FFE::one()])
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![FFE::one(), FFE::one(), FFE::zero()])
        );
        assert_eq!(
            three_bit_iterator.next(),
            Some(vec![FFE::one(), FFE::one(), FFE::one()])
        );
        assert_eq!(three_bit_iterator.next(), None);
    }

    #[test]
    fn test_lagrange_basis_boolean_hypercube_iteration() {
        let minus_one = -FFE::one();
        let check_zero_poly = UniPoly::<FFE>::new(vec![FFE::one(), minus_one]);
        let check_one_poly = UniPoly::<FFE>::new(vec![FFE::zero(), FFE::one()]);
        let mut two_bit_iterator = LagrangeBasisBooleanHyperCube::<FFE>::new(2);
        assert_eq!(two_bit_iterator.total_points, 4);
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
}
