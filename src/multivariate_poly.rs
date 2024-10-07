use std::{
    cmp::Ordering,
    collections::BTreeSet,
    fmt::Debug,
    ops::{Add, Mul, Neg, Sub},
};

use serde::Serialize;

use crate::ff::FiniteFieldElement;

#[derive(Default, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct Var {
    var_index: usize,
    power: usize,
}

impl Var {
    pub fn var_index(&self) -> usize {
        self.var_index
    }

    pub fn power(&self) -> usize {
        self.power
    }
}

impl Debug for Var {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.power == 1 {
            write!(f, "x_{}", self.var_index)?;
        } else {
            write!(f, "(x_{})^{}", self.var_index, self.power)?;
        }
        Ok(())
    }
}

#[derive(Default, Clone, PartialEq, Eq, Serialize)]
pub struct Term<F> {
    coefficient: F,
    vars: Vec<Var>, // x^2y
}

impl<F: FiniteFieldElement> Debug for Term<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.coefficient.is_zero() {
            write!(f, "0")?;
        } else if self.coefficient.is_one() {
            for var in self.vars.iter() {
                if var.power == 1 {
                    write!(f, "x_{}", var.var_index)?;
                } else {
                    write!(f, "(x_{})^{}", var.var_index, var.power)?;
                }
            }
        } else {
            write!(f, "{:?}", self.coefficient)?;
            for var in self.vars.iter() {
                if var.power == 1 {
                    write!(f, "x_{}", var.var_index)?;
                } else {
                    write!(f, "(x_{})^{}", var.var_index, var.power)?;
                }
            }
        }
        Ok(())
    }
}

impl<F: FiniteFieldElement> Term<F> {
    fn is_zero(&self) -> bool {
        self.coefficient.is_zero()
    }

    fn is_one(&self) -> bool {
        self.coefficient.is_one() && self.vars.is_empty()
    }

    fn is_constant(&self) -> bool {
        self.vars.is_empty()
    }

    fn zero() -> Self {
        Term {
            coefficient: F::zero(),
            vars: vec![],
        }
    }

    fn degree(&self) -> usize {
        self.vars.iter().fold(0, |sum, var| sum + var.power)
    }

    fn sort(&mut self) {
        self.vars.sort_by(|a, b| a.var_index.cmp(&b.var_index));
    }

    pub fn coefficient(&self) -> &F {
        &self.coefficient
    }

    pub fn vars(&self) -> &Vec<Var> {
        &self.vars
    }
}

impl<F: FiniteFieldElement> PartialOrd for Term<F> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.degree() != other.degree() {
            return Some(self.degree().cmp(&other.degree()));
        } else {
            for (cur, other) in self.vars.iter().zip(other.vars.iter()) {
                if cur.var_index == other.var_index {
                    if cur.power != other.power {
                        return Some(cur.power.cmp(&other.power));
                    }
                } else {
                    return Some(other.var_index.cmp(&other.var_index));
                }
            }
            Some(Ordering::Equal)
        }
    }
}

impl<F: FiniteFieldElement> Ord for Term<F> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

#[derive(Default, Clone, PartialEq, Eq, Serialize)]
pub struct MultivariatePoly<F> {
    terms: Vec<Term<F>>,
}

impl<F: FiniteFieldElement> Debug for MultivariatePoly<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, term) in self.terms.iter().enumerate() {
            if i == 0 {
                write!(f, "{:?}", term)?;
            } else {
                write!(f, " + {:?}", term)?;
            }
        }
        Ok(())
    }
}

pub trait Polynomial<F>: Sized {
    fn is_zero(&self) -> bool;

    fn is_one(&self) -> bool;

    fn is_constant(&self) -> bool;

    fn new(terms: Vec<Term<F>>) -> Self;

    fn zero() -> Self;

    fn one() -> Self;

    fn terms(&self) -> &Vec<Term<F>>;
}

pub trait MultivariatePolynomial<F>: Polynomial<F> {
    fn remove_zeros(&mut self);
    fn combine_terms(&mut self);
    fn scalar_mul(&mut self, scalar: F);
    fn interpolate(group_of_evaluation_points: &Vec<Vec<F>>, y_values: &Vec<F>) -> Self;
    fn get_lagrange_polynomial_for_single(
        index: usize,
        current_evaluation_point: &F,
        evaluation_points: &BTreeSet<F>,
    ) -> Self;
    fn get_lagrange_polynomial_for_group(
        unique_evaluation_points: &BTreeSet<F>,
        group_of_evaluation_points: &Vec<F>,
    ) -> Self;
}

impl<F: FiniteFieldElement> Polynomial<F> for MultivariatePoly<F> {
    fn new(terms: Vec<Term<F>>) -> Self {
        MultivariatePoly { terms }
    }

    fn is_zero(&self) -> bool {
        self.terms.len() == 1 && self.terms()[0].is_zero()
    }

    fn is_one(&self) -> bool {
        self.terms.len() == 1 && self.terms()[0].is_one()
    }

    fn is_constant(&self) -> bool {
        self.terms.len() == 1 && self.terms()[0].is_constant()
    }

    fn one() -> Self {
        let term = Term {
            coefficient: F::one(),
            vars: vec![],
        };
        let terms = vec![term];
        MultivariatePoly { terms }
    }

    fn zero() -> Self {
        let term = Term {
            coefficient: F::zero(),
            vars: vec![],
        };
        let terms = vec![term];
        MultivariatePoly { terms }
    }

    fn terms(&self) -> &Vec<Term<F>> {
        &self.terms
    }
}

impl<F: FiniteFieldElement + Clone + Neg<Output = F> + Sub<Output = F> + Add<Output = F>>
    MultivariatePolynomial<F> for MultivariatePoly<F>
{
    fn combine_terms(&mut self) {
        let mut terms_dedup: Vec<Term<F>> = Vec::new();

        for term in self.terms.iter() {
            if let Some(prev) = terms_dedup.last_mut() {
                if prev.vars == term.vars {
                    prev.coefficient += term.coefficient.clone();
                    continue;
                }
            }
            terms_dedup.push(term.clone());
        }

        self.terms = terms_dedup
    }

    fn remove_zeros(&mut self) {
        self.terms.retain(|term| !term.coefficient.is_zero());
    }

    fn scalar_mul(&mut self, scalar: F) {
        self.terms = self
            .terms
            .iter()
            .map(|term| Term {
                coefficient: term.coefficient.clone() * scalar.clone(),
                ..term.clone()
            })
            .collect::<Vec<_>>();
    }

    fn interpolate(
        collection_of_group_of_evaluation_points: &Vec<Vec<F>>,
        y_values: &Vec<F>,
    ) -> Self {
        let mut unique_evaluation_points = BTreeSet::new();
        for group in collection_of_group_of_evaluation_points.iter() {
            for evaluation_point in group {
                unique_evaluation_points.insert(evaluation_point.clone());
            }
        }

        let mut resulting_polynomial = Self::zero();
        for (group_of_evaluation_points, y) in collection_of_group_of_evaluation_points
            .iter()
            .zip(y_values)
        {
            let mut group_lagrange_polynomial = Self::get_lagrange_polynomial_for_group(
                &unique_evaluation_points,
                group_of_evaluation_points,
            );
            group_lagrange_polynomial.scalar_mul(y.clone());
            resulting_polynomial = &resulting_polynomial + &group_lagrange_polynomial;
        }
        resulting_polynomial
    }

    fn get_lagrange_polynomial_for_single(
        index: usize,
        current_evaluation_point: &F,
        evaluation_points: &BTreeSet<F>,
    ) -> Self {
        let mut resulting_polynomial = MultivariatePoly::one();
        let mut running_denominator = F::one();

        for (_, evaluation_point) in evaluation_points.iter().enumerate() {
            if current_evaluation_point == evaluation_point {
                continue;
            }
            let var = Var {
                var_index: index,
                power: 1,
            }; // x
            let term = Term {
                coefficient: F::one(),
                vars: vec![var],
            }; // [1x_{i}]
            let constant_term = Term {
                coefficient: -evaluation_point.clone(),
                vars: vec![],
            }; // [-c]
            let terms = vec![term, constant_term];
            let mut numerator = MultivariatePoly::new(terms);
            let denominator = current_evaluation_point.clone() - evaluation_point.clone();
            running_denominator *= denominator.clone();
            numerator.scalar_mul(denominator.inverse().unwrap());
            resulting_polynomial = &resulting_polynomial * &numerator;
        }
        resulting_polynomial
    }

    fn get_lagrange_polynomial_for_group(
        unique_evaluation_points: &BTreeSet<F>,
        group_of_evaluation_points: &Vec<F>,
    ) -> Self {
        let mut group_lagrange_polynomial = Self::one();
        for (index, evaluation) in group_of_evaluation_points.iter().enumerate() {
            let lagrange_polynomial = Self::get_lagrange_polynomial_for_single(
                index + 1,
                evaluation,
                unique_evaluation_points,
            );
            group_lagrange_polynomial = &group_lagrange_polynomial * &lagrange_polynomial;
        }
        group_lagrange_polynomial
    }
}

pub enum ResultOfVarMul {
    Single(Var),
    Vector(Vec<Var>),
}

impl Mul for &Var {
    type Output = ResultOfVarMul;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.var_index == rhs.var_index {
            ResultOfVarMul::Single(Var {
                var_index: self.var_index,
                power: self.power + rhs.power,
            })
        } else {
            let vars = vec![*self, *rhs];
            ResultOfVarMul::Vector(vars)
        }
    }
}

impl<F: FiniteFieldElement + Clone> Mul for &Term<F> {
    type Output = Term<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.is_zero() || rhs.is_zero() {
            return Term::zero();
        } else if self.is_one() {
            return rhs.clone();
        } else if rhs.is_one() {
            return self.clone();
        } else if self.is_constant() {
            return Term {
                coefficient: self.coefficient.clone() * rhs.coefficient.clone(),
                ..rhs.clone()
            };
        } else if rhs.is_constant() {
            return Term {
                coefficient: self.coefficient.clone() * rhs.coefficient.clone(),
                ..self.clone()
            };
        } else {
            let mut resulting_vars = vec![];
            for var_i in self.vars.iter() {
                for var_j in rhs.vars.iter() {
                    let var_mul = var_i * var_j;
                    match var_mul {
                        ResultOfVarMul::Single(var) => {
                            resulting_vars.push(var);
                        }
                        ResultOfVarMul::Vector(vars) => {
                            resulting_vars.extend(vars);
                        }
                    }
                }
            }
            let coefficient = self.coefficient.clone() * rhs.coefficient.clone();
            let mut term = Term {
                coefficient,
                vars: resulting_vars,
            };
            term.sort();
            term
        }
    }
}

pub enum ResultOfTermAdd<F> {
    Single(Term<F>),
    Vector(Vec<Term<F>>),
}

impl<F: FiniteFieldElement + Clone + Add<Output = F>> Add for &Term<F> {
    type Output = ResultOfTermAdd<F>;

    fn add(self, rhs: Self) -> Self::Output {
        if self.vars == rhs.vars {
            ResultOfTermAdd::Single(Term {
                coefficient: self.coefficient.clone() + rhs.coefficient.clone(),
                ..self.clone()
            })
        } else {
            let mut terms = vec![];
            terms.push(self.clone());
            terms.push(rhs.clone());
            ResultOfTermAdd::Vector(terms)
        }
    }
}

impl<F: FiniteFieldElement + Clone + Sub<Output = F> + Neg<Output = F> + Add<Output = F>> Mul
    for &MultivariatePoly<F>
{
    type Output = MultivariatePoly<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.is_zero() || rhs.is_zero() {
            return MultivariatePoly::zero();
        }
        if self.is_one() {
            return rhs.clone();
        }
        if rhs.is_one() {
            return self.clone();
        }
        let mut resulting_terms = vec![];
        for term_i in self.terms() {
            for term_j in rhs.terms() {
                let term_mul = term_i * term_j;
                resulting_terms.push(term_mul);
            }
        }
        let mut poly = MultivariatePoly {
            terms: resulting_terms,
        };
        poly.remove_zeros();
        poly.combine_terms();
        poly
    }
}

impl<F: FiniteFieldElement + Clone + Add<Output = F>> Add for &MultivariatePoly<F> {
    type Output = MultivariatePoly<F>;

    fn add(self, rhs: Self) -> Self::Output {
        if self.is_zero() {
            return rhs.clone();
        }
        if rhs.is_zero() {
            return self.clone();
        }
        let mut result = Vec::new();
        let mut cur_iter = self.terms.iter().peekable();
        let mut other_iter = rhs.terms.iter().peekable();

        loop {
            let which = match (cur_iter.peek(), other_iter.peek()) {
                (Some(cur), Some(other)) => Some(cur.cmp(other)),
                (Some(_), None) => Some(Ordering::Less),
                (None, Some(_)) => Some(Ordering::Greater),
                (None, None) => None,
            };

            let smallest = match which {
                Some(Ordering::Less) => cur_iter.next().unwrap().clone(),
                Some(Ordering::Equal) => {
                    let other = other_iter.next().unwrap();
                    let cur = cur_iter.next().unwrap();
                    Term {
                        coefficient: cur.coefficient.clone() + other.coefficient.clone(),
                        vars: cur.vars.clone(),
                    }
                }
                Some(Ordering::Greater) => other_iter.next().unwrap().clone(),
                None => break,
            };

            result.push(smallest);
        }
        result.retain(|c| !c.is_zero());
        MultivariatePoly::new(result)
    }
}

mod tests {
    use num_bigint::BigInt;

    use crate::ff::{FiniteFieldElement, FFE};

    use super::{MultivariatePoly, MultivariatePolynomial, Polynomial, Term, Var};

    #[test]
    fn test_add() {
        let modulus = BigInt::from(17);
        // x + y + z
        let x = Term {
            coefficient: FFE::new(&BigInt::from(1), &modulus),
            vars: vec![Var {
                var_index: 1,
                power: 1,
            }],
        };
        let y = Term {
            coefficient: FFE::new(&BigInt::from(1), &modulus),
            vars: vec![Var {
                var_index: 2,
                power: 1,
            }],
        };
        let z = Term {
            coefficient: FFE::new(&BigInt::from(1), &modulus),
            vars: vec![Var {
                var_index: 3,
                power: 1,
            }],
        };
        let terms = vec![x.clone(), y.clone(), z.clone()];
        let poly_1 = MultivariatePoly::new(terms);

        // 3x^2 + 5y + 8z^2
        let three_x_square = Term {
            coefficient: FFE::new(&BigInt::from(3), &modulus),
            vars: vec![Var {
                var_index: 1,
                power: 2,
            }],
        };
        let five_y = Term {
            coefficient: FFE::new(&BigInt::from(5), &modulus),
            vars: vec![Var {
                var_index: 2,
                power: 1,
            }],
        };
        let eight_z_square = Term {
            coefficient: FFE::new(&BigInt::from(8), &modulus),
            vars: vec![Var {
                var_index: 3,
                power: 2,
            }],
        };
        let terms = vec![
            three_x_square.clone(),
            five_y.clone(),
            eight_z_square.clone(),
        ];
        let poly_2 = MultivariatePoly::new(terms);

        let actual_poly = &poly_1 + &poly_2;

        let expected_terms = vec![x, y, z, five_y, three_x_square, eight_z_square];
        let expected_poly = MultivariatePoly::new(expected_terms);

        // assert_eq!(actual_poly, expected_poly)
    }

    #[test]
    fn test_interpolate() {
        let modulus = BigInt::from(17);
        let evaluations_points = vec![
            vec![
                FFE::new(&BigInt::from(2), &modulus),
                FFE::new(&BigInt::from(3), &modulus),
            ],
            vec![
                FFE::new(&BigInt::from(2), &modulus),
                FFE::new(&BigInt::from(5), &modulus),
            ],
            vec![
                FFE::new(&BigInt::from(4), &modulus),
                FFE::new(&BigInt::from(8), &modulus),
            ],
        ];
        let evaluations = vec![
            FFE::new(&BigInt::from(5), &modulus),
            FFE::new(&BigInt::from(10), &modulus),
            FFE::new(&BigInt::from(99), &modulus),
        ];
        let poly = MultivariatePoly::interpolate(&evaluations_points, &evaluations);
    }
}
