use std::ops::{Add, Mul, Neg, Sub};

use num_bigint::BigInt;
use num_traits::{One, Zero};

use crate::ff::FiniteFieldElement;

#[derive(Debug, Default, Clone, Copy, PartialEq)]
struct Var {
    var_index: usize,
    power: usize,
}

#[derive(Debug, Default, Clone)]
pub struct Term<F> {
    coefficient: F,
    vars: Vec<Var>, // x^2y
}

impl<F: FiniteFieldElement> Term<F> {
    fn new(coefficient: F, vars: &Vec<Var>) -> Self {
        Term {
            coefficient,
            vars: vars.clone(),
        }
    }

    fn is_zero(&self) -> bool {
        self.coefficient.is_zero() && self.vars.is_empty()
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
}

#[derive(Debug, Default, Clone)]
pub struct MultivariatePoly<F> {
    terms: Vec<Term<F>>,
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
    fn scalar_mul(&mut self, scalar: F);
    fn interpolate(
        num_of_vars: usize,
        group_of_evaluation_points: &Vec<Vec<F>>,
        y_values: &Vec<F>,
    ) -> Result<Self, String>;
    fn get_lagrange_polynomial_for_single(
        current_evaluation_point: &F,
        evaluation_points: &Vec<F>,
    ) -> Self;
    fn get_lagrange_polynomial_for_group(group_of_evaluation_points: &Vec<Self>) -> Self;
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

impl<F: FiniteFieldElement + Clone + Neg<Output = F> + Sub<Output = F>> MultivariatePolynomial<F>
    for MultivariatePoly<F>
{
    fn scalar_mul(&mut self, scalar: F) {
        let mut terms = vec![];
        for term in self.terms.iter() {
            let new_coefficient = term.coefficient.clone() * scalar.clone();
            let new_term = Term::new(new_coefficient, &term.vars);
            terms.push(new_term);
        }
        let new_poly = MultivariatePoly::new(terms);
        *self = new_poly;
    }

    fn interpolate(
        num_of_vars: usize,
        group_of_evaluation_points: &Vec<Vec<F>>,
        y_values: &Vec<F>,
    ) -> Result<Self, String> {
        for group in group_of_evaluation_points {
            if group.len() != num_of_vars {
                return Err(String::from("Invalid collection of vars"));
            }
        }
        let mut checkers_groups = vec![];
        for group in group_of_evaluation_points {
            let mut checkers = vec![];
            for (i, evaluation_point) in group.iter().enumerate() {
                let var = Var {
                    var_index: i + 1,
                    power: 1,
                }; // x_{i}
                let term = Term {
                    coefficient: F::one(),
                    vars: vec![var],
                }; // [1x_{i}]
                let constant_term = Term {
                    coefficient: evaluation_point.clone(),
                    vars: vec![],
                }; // [c]
                let terms = vec![term, constant_term];
                let checker = MultivariatePoly { terms };
                checkers.push(checker);
            }
            checkers_groups.push(checkers);
        }
        let mut resulting_poly = MultivariatePoly::zero();
        for (group, y) in checkers_groups.iter().zip(y_values) {
            let mut one_poly = MultivariatePoly::one();
            for poly in group.iter() {
                one_poly = &one_poly * poly;
            }
            one_poly.scalar_mul(y.clone());
            resulting_poly = &resulting_poly + &one_poly;
        }
        Ok(resulting_poly)
    }

    fn get_lagrange_polynomial_for_single(
        current_evaluation_point: &F,
        evaluation_points: &Vec<F>,
    ) -> Self {
        let mut resulting_polynomial = MultivariatePoly::one();
        let mut running_denominator = F::one();
        for (i, evaluation_point) in evaluation_points.iter().enumerate() {
            if current_evaluation_point == evaluation_point {
                continue;
            }
            let var = Var {
                var_index: i + 1,
                power: 1,
            }; // x_{i}
            let term = Term {
                coefficient: -F::one(),
                vars: vec![var],
            }; // [-1x_{i}]
            let constant_term = Term {
                coefficient: evaluation_point.clone(),
                vars: vec![],
            }; // [c]
            let terms = vec![term, constant_term];
            let numerator = MultivariatePoly::new(terms);
            resulting_polynomial = &resulting_polynomial * &numerator;
            let denominator = current_evaluation_point.clone() - evaluation_point.clone();
            running_denominator *= denominator;
        }
        todo!()
    }

    fn get_lagrange_polynomial_for_group(group_of_evaluation_points: &Vec<Self>) -> Self {
        todo!()
    }
}

enum ResultOfVarMul {
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
            Term {
                coefficient,
                vars: resulting_vars,
            }
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

impl<F: FiniteFieldElement + Clone> Mul for &MultivariatePoly<F> {
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
        MultivariatePoly {
            terms: resulting_terms,
        }
    }
}

impl<F: FiniteFieldElement + Clone> Add for &MultivariatePoly<F> {
    type Output = MultivariatePoly<F>;

    fn add(self, rhs: Self) -> Self::Output {
        if self.is_zero() {
            return rhs.clone();
        }
        if rhs.is_zero() {
            return self.clone();
        }
        let mut self_terms = self.terms.clone();
        let rhs_terms = rhs.terms.clone();
        self_terms.extend(rhs_terms);
        MultivariatePoly::new(self_terms)
        // TODO: Update to solve cases where
        // they might be terms that have same variable combination
    }
}

mod tests {
    use num_bigint::BigInt;

    use crate::ff::{FiniteFieldElement, FFE};

    use super::{MultivariatePoly, MultivariatePolynomial};

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
                FFE::new(&BigInt::from(6), &modulus),
                FFE::new(&BigInt::from(7), &modulus),
            ],
        ];
        let evaluations = vec![
            FFE::new(&BigInt::from(5), &modulus),
            FFE::new(&BigInt::from(5), &modulus),
        ];
        let poly = MultivariatePoly::interpolate(2, &evaluations_points, &evaluations);
        assert!(poly.is_err());

        let evaluations_points = vec![
            vec![
                FFE::new(&BigInt::from(2), &modulus),
                FFE::new(&BigInt::from(3), &modulus),
            ],
            vec![
                FFE::new(&BigInt::from(2), &modulus),
                FFE::new(&BigInt::from(6), &modulus),
            ],
        ];
        let evaluations = vec![
            FFE::new(&BigInt::from(5), &modulus),
            FFE::new(&BigInt::from(5), &modulus),
        ];
        let poly = MultivariatePoly::interpolate(3, &evaluations_points, &evaluations);
        assert!(poly.is_err());

        let evaluations_points = vec![
            vec![
                FFE::new(&BigInt::from(2), &modulus),
                FFE::new(&BigInt::from(3), &modulus),
            ],
            vec![
                FFE::new(&BigInt::from(2), &modulus),
                FFE::new(&BigInt::from(5), &modulus),
            ],
        ];
        let evaluations = vec![
            FFE::new(&BigInt::from(5), &modulus),
            FFE::new(&BigInt::from(10), &modulus),
        ];
        let poly = MultivariatePoly::interpolate(2, &evaluations_points, &evaluations);
        assert!(poly.is_ok());
        println!("{:?}", poly.unwrap());
    }
}
