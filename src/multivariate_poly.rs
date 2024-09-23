use std::{
    collections::BTreeSet,
    ops::{Add, Mul},
};

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

    fn one() -> Self {
        Term {
            coefficient: F::one(),
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
        field: BigInt,
    ) -> Result<Self, String>;
}

impl<F: FiniteFieldElement> Polynomial<F> for MultivariatePoly<F> {
    fn is_zero(&self) -> bool {
        self.terms.len() == 1
            && self.terms[0].coefficient.element() == BigInt::zero()
            && self.terms[0].vars.is_empty()
    }

    fn is_one(&self) -> bool {
        self.terms.len() == 1
            && self.terms[0].coefficient.element() == BigInt::one()
            && self.terms[0].vars.is_empty()
    }

    fn new(terms: Vec<Term<F>>) -> Self {
        MultivariatePoly { terms }
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

impl<F: FiniteFieldElement + Clone> MultivariatePolynomial<F> for MultivariatePoly<F> {
    fn scalar_mul(&mut self, scalar: F) {
        todo!()
    }

    fn interpolate(
        num_of_vars: usize,
        group_of_evaluation_points: &Vec<Vec<F>>,
        y_values: &Vec<F>,
        field: BigInt,
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
        todo!()
    }
}

mod tests {}
