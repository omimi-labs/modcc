use std::{
    cmp::Ordering,
    collections::BTreeSet,
    fmt::Debug,
    ops::{Add, Mul, Neg, Sub},
};

use num_bigint::BigInt;
use num_traits::One;
use serde::Serialize;

use crate::ff::FiniteFieldElement;

#[derive(Default, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct Var {
    index: usize,
    power: usize,
}

impl Var {
    pub fn new(index: usize, power: usize) -> Self {
        Var { index, power }
    }

    pub fn index(&self) -> usize {
        self.index
    }

    pub fn power(&self) -> usize {
        self.power
    }

    pub fn pow(&self, n: usize) -> Self {
        Var::new(self.index, self.power * n)
    }

    fn to_latex(&self) -> String {
        let latex: String;
        if self.power() == 1 {
            latex = format!("x_{{{}}}", self.index())
        } else {
            latex = format!("x_{{{}}}^{{{}}}", self.index(), self.power())
        }
        latex
    }

    fn from_latex(latex_string: String) -> Self {
        let power_and_index: Vec<&str> = latex_string.split(&['{', '}', '^'][..]).collect();
        if power_and_index.len() == 3 {
            let index: usize = power_and_index[1].parse().unwrap();
            return Self::new(index, 1);
        } else {
            let index: usize = power_and_index[1].parse().unwrap();
            let power: usize = power_and_index[4].parse().unwrap();
            return Self::new(index, power);
        }
    }
}

impl Debug for Var {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.power == 1 {
            write!(f, "x_{}", self.index)?;
        } else {
            write!(f, "(x_{})^{}", self.index, self.power)?;
        }
        Ok(())
    }
}

#[derive(Default, Clone, PartialEq, Eq, Serialize)]
pub struct Term<F> {
    coefficient: F,
    vars: Vec<Var>, // x^2y
}

impl<F: FiniteFieldElement + Clone + Add<Output = F>> Debug for Term<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_zero() {
            write!(f, "0")?;
        } else if self.is_one() {
            if self.is_one_and_constant() {
                write!(f, "1")?;
            } else {
                for var in self.vars.iter() {
                    write!(f, "{:?}", var)?;
                }
            }
        } else if self.is_constant() {
            write!(f, "{:?}", self.coefficient)?;
        } else {
            write!(f, "{:?}", self.coefficient)?;
            for var in self.vars.iter() {
                write!(f, "{:?}", var)?;
            }
        }
        Ok(())
    }
}

impl<F: FiniteFieldElement + Clone + Add<Output = F>> Term<F> {
    pub fn new(coefficient: F, vars: Vec<Var>) -> Self {
        Term { coefficient, vars }
    }

    fn is_zero(&self) -> bool {
        self.coefficient.is_zero()
    }

    fn is_one(&self) -> bool {
        self.coefficient.is_one()
    }

    fn is_constant(&self) -> bool {
        self.vars.is_empty()
    }

    fn is_one_and_constant(&self) -> bool {
        self.is_one() && self.is_constant()
    }

    pub fn one() -> Self {
        Term {
            coefficient: F::one(),
            vars: vec![],
        }
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
        self.vars.sort_by(|a, b| a.power.cmp(&b.power));
        self.vars.sort_by(|a, b| a.index.cmp(&b.index));
    }

    pub fn coefficient(&self) -> &F {
        &self.coefficient
    }

    pub fn vars(&self) -> &Vec<Var> {
        &self.vars
    }

    fn to_latex(&self) -> String {
        let mut latex: String;
        if self.is_zero() {
            latex = String::from("0");
            return latex;
        } else if self.is_one_and_constant() {
            latex = String::from("1");
            return latex;
        } else if self.is_one() {
            latex = String::from("");
            for var in self.vars() {
                let var_latex = var.to_latex();
                latex += &var_latex;
            }
            return latex;
        } else if self.is_constant() {
            latex = self.coefficient().element().to_string();
            return latex;
        } else {
            latex = self.coefficient().element().to_string();
            for var in self.vars() {
                let var_latex = var.to_latex();
                latex += &var_latex;
            }
            return latex;
        }
    }

    fn from_latex(latex_string: String, modulus: &BigInt) -> Self {
        let chars: Vec<&str> = latex_string.split("x").collect();
        if chars.len() == 1 {
            let constant: BigInt = chars[0].parse().unwrap();
            let coefficient = F::new(&constant, modulus);
            let term = Self::new(coefficient, vec![]);
            return term;
        } else if chars[0] == "" {
            let constant: BigInt = BigInt::one();
            let coefficient = F::new(&constant, modulus);
            let mut vars = vec![];
            for (i, var_str) in chars.iter().enumerate() {
                if i == 0 {
                    continue;
                }
                let var_latex = format!("x{}", var_str);
                let var = Var::from_latex(var_latex);
                vars.push(var);
            }
            let term = Self::new(coefficient, vars);
            return term;
        } else {
            let constant: BigInt = chars[0].parse().unwrap();
            let coefficient = F::new(&constant, modulus);
            let mut vars = vec![];
            for (i, var_str) in chars.iter().enumerate() {
                if i == 0 {
                    continue;
                }
                let var_latex = format!("x{}", var_str);
                let var = Var::from_latex(var_latex);
                vars.push(var);
            }
            let term = Self::new(coefficient, vars);
            return term;
        }
    }

    // Sums the powers the duplicated variables.
    pub fn combine_terms(&mut self) {
        self.sort();
        let mut vars_dedup: Vec<Var> = Vec::new();

        for var in self.vars().iter() {
            if let Some(prev) = vars_dedup.last_mut() {
                if prev.index() == var.index() {
                    prev.power += var.power();
                    continue;
                }
            }
            vars_dedup.push(*var);
        }
        self.vars = vars_dedup
    }

    pub fn evaluate(&self, evaluation_points: &Vec<(usize, F)>) -> Self {
        let mut term = self.clone();
        let vars = term.vars().clone();
        for var in vars.iter() {
            if evaluation_points
                .iter()
                .any(|(index, _)| index == &var.index())
            {
                let index = evaluation_points
                    .iter()
                    .position(|(i, _)| i == &var.index());
                let evaluation =
                    (evaluation_points[index.unwrap()].1).pow(var.power.try_into().unwrap());
                term.coefficient *= evaluation;
                let index = term.vars().iter().position(|v| v.index() == var.index());
                term.vars.remove(index.unwrap());
            }
        }
        term.combine_terms();
        term
    }
}

impl<F: FiniteFieldElement + Clone + Add<Output = F>> PartialOrd for Term<F> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.degree() != other.degree() {
            return Some(self.degree().cmp(&other.degree()));
        } else {
            for (cur, other) in self.vars.iter().zip(other.vars.iter()) {
                if cur.index == other.index {
                    if cur.power != other.power {
                        return Some(cur.power.cmp(&other.power));
                    }
                } else {
                    return Some(cur.index.cmp(&other.index));
                }
            }
            Some(Ordering::Equal)
        }
    }
}

impl<F: FiniteFieldElement + Clone + Add<Output = F>> Ord for Term<F> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl<F: FiniteFieldElement + Clone + Add<Output = F>> Mul for &Term<F> {
    type Output = Term<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.is_zero() || rhs.is_zero() {
            return Term::zero();
        } else if self.is_one_and_constant() {
            return rhs.clone();
        } else if rhs.is_one_and_constant() {
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
            let mut self_vars = self.vars().clone();
            self_vars.extend(rhs.vars());
            let coefficient = self.coefficient.clone() * rhs.coefficient.clone();
            let mut term = Term {
                coefficient,
                vars: self_vars,
            };
            term.combine_terms();
            term
        }
    }
}

#[derive(Debug, Serialize)]
pub struct InterpolationStep {
    step_1: String,
    step_2: String,
    step_3: String,
    step_4: String,
}

#[derive(Debug, Serialize)]
pub struct Steps {
    interpolation_steps: Vec<InterpolationStep>,
    combination_step: String,
}

#[derive(Debug, Serialize)]
pub struct Properties {
    num_of_vars: usize,
    num_of_evaluations: usize,
}

#[derive(Default, Clone, PartialEq, Eq, Serialize)]
pub struct MultivariatePoly<F> {
    terms: Vec<Term<F>>,
}

impl<F: FiniteFieldElement + Clone + Add<Output = F>> Debug for MultivariatePoly<F> {
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
    fn create_latex_lhs(num_of_vars: usize) -> String;
    fn create_latex_rhs(&self) -> String;
    fn to_latex(&self, num_of_vars: usize) -> String;
    fn from_latex(latex_string: String, modulus: &BigInt) -> Self;
    fn sort(&mut self);
    fn remove_zeros(&mut self);
    fn combine_terms(&mut self);
    fn scalar_mul(&mut self, scalar: F);
    fn interpolate(
        group_of_evaluation_points: &Vec<Vec<F>>,
        y_values: &Vec<F>,
    ) -> (Self, Steps, Properties);
    fn get_lagrange_polynomial_for_single(
        index: usize,
        current_evaluation_point: &F,
        evaluation_points: &BTreeSet<F>,
    ) -> Self;
    fn get_lagrange_polynomial_for_group(
        unique_evaluation_points: &BTreeSet<F>,
        group_of_evaluation_points: &Vec<F>,
    ) -> (Self, InterpolationStep);
    fn evaluation(&self, evaluation_points: &Vec<(usize, F)>) -> Self;
    fn full_evaluation(&self, evaluation_points: &Vec<(usize, F)>) -> F;
    fn partial_evaluation(&self, evaluation_points: &Vec<(usize, F)>) -> Self;
}

impl<F: FiniteFieldElement + Clone + Add<Output = F>> Polynomial<F> for MultivariatePoly<F> {
    fn new(terms: Vec<Term<F>>) -> Self {
        MultivariatePoly { terms }
    }

    fn is_zero(&self) -> bool {
        self.terms.len() == 1 && self.terms()[0].is_zero()
    }

    fn is_one(&self) -> bool {
        self.terms.len() == 1 && self.terms()[0].is_one_and_constant()
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
    fn create_latex_lhs(num_of_vars: usize) -> String {
        let mut value = String::from("f(");
        for i in 0..num_of_vars {
            if i == num_of_vars - 1 {
                value += &format!("x_{}", i + 1);
            } else {
                value += &format!("x_{},", i + 1);
            }
        }
        value += ")";
        value
    }

    fn create_latex_rhs(&self) -> String {
        let mut latex: String;
        if self.is_zero() || self.is_one() || self.is_constant() {
            latex = self.terms[0].to_latex();
            return latex;
        } else if self.terms().len() == 1 {
            latex = self.terms[0].to_latex();
            return latex;
        } else {
            if self.terms[0].is_constant() {
                let constant = self.terms[0].to_latex();
                latex = self.terms[1].to_latex();
                for (i, term) in self.terms().iter().enumerate() {
                    if i == 0 || i == 1 {
                        continue;
                    }
                    let term_latex = term.to_latex();
                    latex += &format!(" + {}", term_latex);
                }
                latex += &format!(" + {}", constant);
                return latex;
            } else {
                latex = self.terms[0].to_latex();
                for (i, term) in self.terms().iter().enumerate() {
                    if i == 0 {
                        continue;
                    }
                    let term_latex = term.to_latex();
                    latex += &format!(" + {}", term_latex);
                }
                return latex;
            }
        }
    }

    fn to_latex(&self, num_of_vars: usize) -> String {
        format!(
            "${} = {}$",
            Self::create_latex_lhs(num_of_vars),
            self.create_latex_rhs(),
        )
    }

    fn from_latex(latex_string: String, modulus: &BigInt) -> Self {
        let rhs: Vec<&str> = latex_string.split("=").collect();
        let terms_strs: Vec<&str> = rhs[1].split("+").collect();
        let trimmed_term_strs: Vec<&str> =
            terms_strs.iter().map(|term_str| term_str.trim()).collect();
        let mut terms = vec![];
        for term_str in trimmed_term_strs {
            let term = Term::<F>::from_latex(term_str.to_string(), modulus);
            terms.push(term);
        }
        if terms.len() != 1 && terms.last().unwrap().is_constant() {
            let last_term = terms.pop().unwrap();
            terms.insert(0, last_term);
            let poly = Self::new(terms);
            return poly;
        } else {
            let poly = Self::new(terms);
            return poly;
        }
    }

    fn sort(&mut self) {
        self.terms.sort_by(|a, b| a.cmp(&b));
    }

    fn remove_zeros(&mut self) {
        self.terms.retain(|term| !term.is_zero());
    }

    fn combine_terms(&mut self) {
        let mut terms_dedup: Vec<Term<F>> = Vec::new();

        for term in self.terms.iter() {
            if let Some(prev) = terms_dedup.last_mut() {
                if prev.vars() == term.vars() {
                    prev.coefficient += term.coefficient.clone();
                    continue;
                }
            }
            terms_dedup.push(term.clone());
        }

        self.terms = terms_dedup
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

    fn evaluation(&self, evaluation_points: &Vec<(usize, F)>) -> Self {
        let mut poly = self.clone();
        for term in poly.terms.iter_mut() {
            let new_term = term.evaluate(evaluation_points);
            *term = new_term;
        }
        poly.sort();
        poly.combine_terms();
        poly.remove_zeros();
        poly
    }

    fn full_evaluation(&self, evaluation_points: &Vec<(usize, F)>) -> F {
        let poly = self.evaluation(evaluation_points);
        if poly.terms().len() == 0 {
            F::zero()
        } else {
            poly.terms[0].coefficient.clone()
        }
    }

    fn partial_evaluation(&self, evaluation_points: &Vec<(usize, F)>) -> Self {
        let poly = self.evaluation(evaluation_points);
        poly
    }

    /**
     * collection_of_group_of_evaluation_points: ((0, 1), (1, 0), (2, 5), (6, 9))
     * y_values: (1, 3, 4, 5)
     *
     * Steps:
     *      1. Create a set of single evaluation points
     *      2.
     */
    fn interpolate(
        collection_of_group_of_evaluation_points: &Vec<Vec<F>>,
        y_values: &Vec<F>,
    ) -> (Self, Steps, Properties) {
        println!("{:?}", collection_of_group_of_evaluation_points);
        let mut interpolation_steps = vec![];
        let mut combination_step = String::from("");
        let mut unique_evaluation_points = BTreeSet::new();
        for group in collection_of_group_of_evaluation_points.iter() {
            for evaluation_point in group {
                unique_evaluation_points.insert(evaluation_point.clone());
            }
        }
        let mut resulting_polynomial = Self::zero();
        for (i, (group_of_evaluation_points, y)) in collection_of_group_of_evaluation_points
            .iter()
            .zip(y_values)
            .enumerate()
        {
            let (mut group_lagrange_polynomial, step) = Self::get_lagrange_polynomial_for_group(
                &unique_evaluation_points,
                group_of_evaluation_points,
            );
            if i == 0 {
                combination_step +=
                    &format!("{:?}({})", y, group_lagrange_polynomial.create_latex_rhs());
            } else {
                combination_step += &format!(
                    " + {:?}({})",
                    y,
                    group_lagrange_polynomial.create_latex_rhs()
                );
            }
            interpolation_steps.push(step);
            group_lagrange_polynomial.scalar_mul(y.clone());
            resulting_polynomial = &resulting_polynomial + &group_lagrange_polynomial;
        }

        let steps = Steps {
            interpolation_steps,
            combination_step,
        };
        let properties = Properties {
            num_of_vars: collection_of_group_of_evaluation_points[0].len(),
            num_of_evaluations: y_values.len(),
        };
        (resulting_polynomial, steps, properties)
    }

    fn get_lagrange_polynomial_for_group(
        unique_evaluation_points: &BTreeSet<F>,
        group_of_evaluation_points: &Vec<F>,
    ) -> (Self, InterpolationStep) {
        let mut step_1 = String::from("L(");
        let mut step_2 = String::from("");
        let mut step_3 = String::from("");

        let mut group_lagrange_polynomial = Self::one();
        for (index, evaluation) in group_of_evaluation_points.iter().enumerate() {
            if index == 0 {
                step_1 += &format!("{:?}", evaluation);
            } else {
                step_1 += &format!(",{:?}", evaluation);
            }
            step_2 += &format!("(L_{:?})", evaluation);
            let lagrange_polynomial = Self::get_lagrange_polynomial_for_single(
                index + 1,
                evaluation,
                unique_evaluation_points,
            );
            step_3 += &format!("({})", lagrange_polynomial.create_latex_rhs());
            group_lagrange_polynomial = &group_lagrange_polynomial * &lagrange_polynomial;
        }
        step_1 += ")";
        let step_4 = format!("{}", group_lagrange_polynomial.clone().create_latex_rhs());
        let interpolation_step = InterpolationStep {
            step_1,
            step_2,
            step_3,
            step_4,
        };
        (group_lagrange_polynomial, interpolation_step)
    }

    fn get_lagrange_polynomial_for_single(
        index: usize,
        current_evaluation_point: &F,
        evaluation_points: &BTreeSet<F>,
    ) -> Self {
        let mut resulting_polynomial = MultivariatePoly::one();
        for evaluation_point in evaluation_points.iter() {
            if current_evaluation_point == evaluation_point {
                continue;
            }
            let var = Var {
                index: index,
                power: 1,
            }; // x_{index}
            let term = Term {
                coefficient: F::one(),
                vars: vec![var],
            }; // [1x_{i}]
            let constant_term = Term {
                coefficient: -evaluation_point.clone(),
                vars: vec![],
            }; // [-c]

            let terms = vec![constant_term, term];
            let mut numerator = MultivariatePoly::new(terms);
            let denominator = current_evaluation_point.clone() - evaluation_point.clone();
            numerator.scalar_mul(denominator.inverse().unwrap());
            resulting_polynomial = &resulting_polynomial * &numerator;
        }
        resulting_polynomial
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
        let mut terms = vec![];
        for self_term in self.terms.iter() {
            for rhs_term in rhs.terms.iter() {
                let term = self_term * rhs_term;
                terms.push(term);
            }
        }
        let mut poly = MultivariatePoly::new(terms);
        poly.remove_zeros();
        poly.sort();
        poly.combine_terms();
        poly
    }
}

impl<F: FiniteFieldElement + Clone + Add<Output = F> + Neg<Output = F> + Sub<Output = F>> Add
    for &MultivariatePoly<F>
{
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

    #[allow(unused_imports)]
    use num_traits::{One, Zero};
    use rand::{thread_rng, Rng};

    use crate::ff::{FiniteFieldElement, FFE};

    #[allow(unused_imports)]
    use super::{MultivariatePoly, MultivariatePolynomial, Polynomial, Term, Var};

    fn generate_point(num_of_points: usize, modulus: &BigInt) -> Vec<FFE> {
        let mut points = vec![];
        let mut rng = thread_rng();
        for _ in 0..num_of_points {
            let value: usize = rng.gen();
            let element = BigInt::from(value);
            let var = FFE::new(&element, modulus);
            points.push(var);
        }
        points
    }

    fn generate_random_evaluation_points(
        evaluation_points: &mut Vec<Vec<FFE>>,
        num_of_vars: usize,
        num_of_points: usize,
        modulus: &BigInt,
    ) -> Vec<Vec<FFE>> {
        for _ in 0..num_of_points {
            let point = generate_point(num_of_vars, modulus);
            if evaluation_points.contains(&point) {
                let point = generate_point(num_of_vars, modulus);
                evaluation_points.push(point);
            } else {
                evaluation_points.push(point);
            }
        }
        evaluation_points.to_vec()
    }

    fn get_evaluation_point(index: usize, evaluation_points: &Vec<Vec<FFE>>) -> Vec<(usize, FFE)> {
        let point: &Vec<FFE> = &evaluation_points[index];
        let mut evaluation_point = vec![];
        for (i, point) in point.iter().enumerate() {
            evaluation_point.push((i + 1, point.clone()));
        }
        evaluation_point
    }

    #[test]
    fn test_combining_terms() {
        let modulus = BigInt::from(17);

        let x_1 = Var::new(1, 1);
        let x_2 = Var::new(2, 1);
        let x_1_square = Var::new(1, 2);
        let x_2_cube = Var::new(2, 3);

        let coefficient = FFE::new(&BigInt::from(3), &modulus);
        let vars = vec![x_1, x_2, x_1_square, x_2_cube];
        let mut term = Term::new(coefficient, vars);
        term.combine_terms();
        assert!(term.vars().len() == 2);
        assert_eq!(term.vars()[0], Var::new(1, 3));
        assert_eq!(term.vars()[1], Var::new(2, 4));

        let x_3_exp_five = Var::new(3, 5);
        let x_3 = Var::new(3, 1);
        let x_4_exp_six = Var::new(4, 6);
        let x_5_exp_2 = Var::new(5, 2);
        let coefficient = FFE::new(&BigInt::from(3), &modulus);
        let vars = vec![x_3_exp_five, x_3, x_4_exp_six, x_5_exp_2];
        let mut term = Term::new(coefficient, vars);
        term.combine_terms();
        assert!(term.vars().len() == 3);
        assert_eq!(term.vars()[0], Var::new(3, 6));
        assert_eq!(term.vars()[1], Var::new(4, 6));
        assert_eq!(term.vars()[2], Var::new(5, 2));

        let coefficient = FFE::new(&BigInt::from(3), &modulus);
        let vars = vec![
            x_1,
            x_3_exp_five,
            x_2,
            x_3,
            x_4_exp_six,
            x_5_exp_2,
            x_1_square,
            x_2_cube,
        ];
        let mut term = Term::new(coefficient, vars);
        term.combine_terms();
        assert!(term.vars().len() == 5);
        assert_eq!(term.vars()[0], Var::new(1, 3));
        assert_eq!(term.vars()[1], Var::new(2, 4));
        assert_eq!(term.vars()[2], Var::new(3, 6));
        assert_eq!(term.vars()[3], Var::new(4, 6));
        assert_eq!(term.vars()[4], Var::new(5, 2));
    }

    #[test]
    fn test_term_comparison() {
        let modulus = BigInt::from(17);

        // Scenario 1: Unequal degrees
        let x_1 = Var::new(1, 1);
        let x_2 = Var::new(2, 1);
        let x_1_square = Var::new(1, 2);
        let x_2_cube = Var::new(2, 3);
        let coefficient = FFE::new(&BigInt::from(3), &modulus);
        let vars = vec![x_1, x_2, x_1_square, x_2_cube];
        let mut term_1 = Term::new(coefficient, vars);
        term_1.combine_terms();

        let x_3_exp_five = Var::new(3, 5);
        let x_3 = Var::new(3, 1);
        let x_4_exp_six = Var::new(4, 6);
        let x_5_exp_2 = Var::new(5, 2);
        let coefficient = FFE::new(&BigInt::from(3), &modulus);
        let vars = vec![x_3_exp_five, x_3, x_4_exp_six, x_5_exp_2];
        let mut term_2 = Term::new(coefficient, vars);
        term_2.combine_terms();

        assert!(term_2 > term_1);

        // Scenario 2.1: Equal degree and Exact variables
        let vars = vec![x_1, x_2, x_1_square, x_2_cube];
        let coefficient = FFE::new(&BigInt::from(3), &modulus);
        let mut term_1 = Term::new(coefficient, vars);
        term_1.combine_terms();

        let vars = vec![x_1, x_2, x_2_cube, x_1_square];
        let coefficient = FFE::new(&BigInt::from(3), &modulus);
        let mut term_2 = Term::new(coefficient, vars);
        term_2.combine_terms();

        assert!(term_1 == term_2);

        // Scenario 2.2: Equal degree and unequal variables
        let vars = vec![x_1, x_2, x_3];
        let coefficient = FFE::new(&BigInt::from(3), &modulus);
        let mut term_1 = Term::new(coefficient, vars);
        term_1.combine_terms();

        let x_4 = Var::new(4, 1);
        let vars = vec![x_1, x_2, x_4];
        let coefficient = FFE::new(&BigInt::from(3), &modulus);
        let mut term_2 = Term::new(coefficient, vars);
        term_2.combine_terms();

        assert!(term_1 < term_2);

        // Scenario 2.3: Equal degree, equal variable and unequal power
        let vars = vec![x_1_square, x_2];
        let coefficient = FFE::new(&BigInt::from(3), &modulus);
        let mut term_1 = Term::new(coefficient, vars);
        term_1.combine_terms();

        let x_2_square = Var::new(2, 2);
        let vars = vec![x_1, x_2_square];
        let coefficient = FFE::new(&BigInt::from(3), &modulus);
        let mut term_2 = Term::new(coefficient, vars);
        term_2.combine_terms();

        assert!(term_1 > term_2);

        // EDGE CASES

        // Unequal number of vars
        let x_1_cube = Var::new(1, 3);
        let vars = vec![x_1_cube];
        let coefficient = FFE::new(&BigInt::from(3), &modulus);
        let term_1 = Term::new(coefficient, vars);

        let vars = vec![x_1_square, x_2];
        let coefficient = FFE::new(&BigInt::from(3), &modulus);
        let term_2 = Term::new(coefficient, vars);

        assert!(term_1 > term_2);
    }

    #[test]
    fn test_term_mul() {
        // x_1
        let x_1 = Var::new(1, 1);
        let term_x_1 = Term::new(FFE::one(), vec![x_1]);
        // x_2
        let x_2 = Var::new(2, 1);
        let term_x_2 = Term::new(FFE::one(), vec![x_2]);

        // x_1x_2
        let term_x_1_x_2 = &term_x_1 * &term_x_2;
        assert_eq!(Term::new(FFE::one(), vec![x_1, x_2]), term_x_1_x_2);

        // (x_1)^2
        let term_x_1_sq = &term_x_1 * &term_x_1;
        assert_eq!(Term::new(FFE::one(), vec![Var::new(1, 2)]), term_x_1_sq);

        // (x_2)^2
        let term_x_2_sq = &term_x_2 * &term_x_2;
        assert_eq!(Term::new(FFE::one(), vec![Var::new(2, 2)]), term_x_2_sq);

        let term_x_1_x_2_sq = &term_x_1 * &term_x_2_sq;
        assert_eq!(
            Term::new(FFE::one(), vec![x_1, Var::new(2, 2)]),
            term_x_1_x_2_sq
        );
        let term_x_2_sq_x_1 = &term_x_2_sq * &term_x_1;
        assert_eq!(
            Term::new(FFE::one(), vec![x_1, Var::new(2, 2)]),
            term_x_2_sq_x_1
        );
        // (x_1)^2(x_2)^2
        let term_x_1_sq_x_2_sq = &term_x_1_sq * &term_x_2_sq;
        assert_eq!(
            Term::new(FFE::one(), vec![Var::new(1, 2), Var::new(2, 2)]),
            term_x_1_sq_x_2_sq
        );
        // (x_1)^2(x_2)^2
        let term_x_2_sq_x_1_sq = &term_x_2_sq * &term_x_1_sq;
        assert_eq!(
            Term::new(FFE::one(), vec![Var::new(1, 2), Var::new(2, 2)]),
            term_x_2_sq_x_1_sq
        );
        // (x_1)^3(x_2)^2
        let term_x_1_cub_x_2_sq = &term_x_1 * &term_x_2_sq_x_1_sq;
        assert_eq!(
            Term::new(FFE::one(), vec![Var::new(1, 3), Var::new(2, 2)]),
            term_x_1_cub_x_2_sq
        );
        // (x_1)^3(x_2)^2
        let term_x_2_sq_x_1_cub = &term_x_2_sq_x_1_sq * &term_x_1;
        assert_eq!(
            Term::new(FFE::one(), vec![Var::new(1, 3), Var::new(2, 2)]),
            term_x_2_sq_x_1_cub
        );
        let modulus = BigInt::from(17);
        let three = FFE::new(&BigInt::from(3), &modulus);
        let constant_term_three = Term::new(three.clone(), vec![]);
        // 3x_1
        let three_x_1 = &constant_term_three * &term_x_1;
        assert_eq!(Term::new(three.clone(), vec![x_1]), three_x_1);

        let x_1_three = &term_x_1 * &constant_term_three;
        assert_eq!(Term::new(three, vec![x_1]), x_1_three);

        // 5x_1x_2^6x_3 * 6x_1x_4x_5^9
        let five = FFE::new(&BigInt::from(5), &modulus);
        let x_2_sixth = x_2.pow(6);
        let x_3 = Var::new(3, 1);
        let x_5 = Var::new(5, 1);
        let vars = vec![x_2_sixth, x_1, x_3, x_5];
        let term_1 = Term::new(five, vars);

        let three = FFE::new(&BigInt::from(3), &modulus);
        let x_4 = Var::new(4, 1);
        let x_5_nineth = x_5.pow(9);
        let vars = vec![x_1, x_4, x_5_nineth];
        let term_2 = Term::new(three, vars);

        let term_3 = &term_1 * &term_2;
        assert_eq!(
            Term::new(
                FFE::new(&BigInt::from(15), &modulus),
                vec![x_1.pow(2), x_2_sixth, x_3, x_4, x_5.pow(10)]
            ),
            term_3
        );
    }

    #[test]
    fn test_poly_mul() {
        // Scenario 1: same poly

        // 5x_1x_2^6x_3 + 3x_1x_4x_5^9
        let modulus = BigInt::from(17);
        // x_1
        let x_1 = Var::new(1, 1);
        // x_2
        let x_2 = Var::new(2, 1);

        let five = FFE::new(&BigInt::from(5), &modulus);
        let x_2_sixth = x_2.pow(6);
        let x_3 = Var::new(3, 1);
        let x_5 = Var::new(5, 1);
        let vars = vec![x_2_sixth, x_1, x_3, x_5];
        let mut term_1 = Term::new(five, vars);
        term_1.combine_terms();

        let three = FFE::new(&BigInt::from(3), &modulus);
        let x_4 = Var::new(4, 1);
        let x_5_nineth = x_5.pow(9);
        let vars = vec![x_1, x_4, x_5_nineth];
        let mut term_2 = Term::new(three, vars);
        term_2.combine_terms();

        let poly_1 = MultivariatePoly::new(vec![term_1, term_2]);

        // [5x_1x_2^6x_3 + 3x_1x_4x_5^9] * [5x_1x_2^6x_3 + 3x_1x_4x_5^9]
        let poly_2 = &poly_1 * &poly_1;
        let terms = vec![
            Term::new(
                FFE::new(&BigInt::from(25), &modulus),
                vec![x_1.pow(2), x_2.pow(12), x_3.pow(2), x_5.pow(2)],
            ),
            Term::new(
                FFE::new(&BigInt::from(30), &modulus),
                vec![x_1.pow(2), x_2.pow(6), x_3, x_4, x_5.pow(10)],
            ),
            Term::new(
                FFE::new(&BigInt::from(9), &modulus),
                vec![x_1.pow(2), x_4.pow(2), x_5.pow(18)],
            ),
        ];
        let expected_poly = MultivariatePoly::new(terms);
        assert_eq!(expected_poly, poly_2);

        // [5x_1x_2^6x_3 + 3x_1x_4x_5^9] * [5x_1x_2^6x_3 + 3x_1x_4x_5^9] * [5x_1x_2^6x_3 + 3x_1x_4x_5^9] * [5x_1x_2^6x_3 + 3x_1x_4x_5^9]
        let poly_3 = &poly_2 * &poly_2;
        let terms = vec![
            Term::new(
                FFE::new(&BigInt::from(625), &modulus),
                vec![x_1.pow(4), x_2.pow(24), x_3.pow(4), x_5.pow(4)],
            ),
            Term::new(
                FFE::new(&BigInt::from(1500), &modulus),
                vec![x_1.pow(4), x_2.pow(18), x_3.pow(3), x_4, x_5.pow(12)],
            ),
            Term::new(
                FFE::new(&BigInt::from(1350), &modulus),
                vec![x_1.pow(4), x_2.pow(12), x_3.pow(2), x_4.pow(2), x_5.pow(20)],
            ),
            Term::new(
                FFE::new(&BigInt::from(540), &modulus),
                vec![x_1.pow(4), x_2.pow(6), x_3, x_4.pow(3), x_5.pow(28)],
            ),
            Term::new(
                FFE::new(&BigInt::from(81), &modulus),
                vec![x_1.pow(4), x_4.pow(4), x_5.pow(36)],
            ),
        ];
        let expected_poly = MultivariatePoly::new(terms);
        assert_eq!(expected_poly, poly_3);

        // Scenario 2: different poly of same number of terms

        // 10x_1x_2x_3
        let ten = FFE::new(&BigInt::from(10), &modulus);
        let term_1 = Term::new(ten, vec![x_1, x_2, x_3]);

        // 11(x_4)^3x_5
        let eleven = FFE::new(&BigInt::from(11), &modulus);
        let term_2 = Term::new(eleven, vec![x_4.pow(3), x_5]);

        let poly_4 = MultivariatePoly::new(vec![term_1, term_2]);

        let poly_5 = &poly_1 * &poly_4;

        let poly_6 = &poly_4 * &poly_1;

        assert_eq!(poly_5, poly_6);

        let terms = vec![
            Term::new(
                FFE::new(&BigInt::from(50), &modulus),
                vec![x_1.pow(2), x_2.pow(7), x_3.pow(2), x_5],
            ),
            Term::new(
                FFE::new(&BigInt::from(55), &modulus),
                vec![x_1, x_2.pow(6), x_3, x_4.pow(3), x_5.pow(2)],
            ),
            Term::new(
                FFE::new(&BigInt::from(30), &modulus),
                vec![x_1.pow(2), x_2, x_3, x_4, x_5.pow(9)],
            ),
            Term::new(
                FFE::new(&BigInt::from(33), &modulus),
                vec![x_1, x_4.pow(4), x_5.pow(10)],
            ),
        ];

        let expected_poly = MultivariatePoly::new(terms);

        assert_eq!(expected_poly, poly_6);

        let poly_7 = &poly_3 * &poly_5;

        let poly_8 = &poly_5 * &poly_3;

        assert_eq!(poly_7, poly_8);

        let terms = vec![
            Term::new(
                FFE::new(&BigInt::from(31250), &modulus),
                vec![x_1.pow(6), x_2.pow(31), x_3.pow(6), x_5.pow(5)],
            ),
            Term::new(
                FFE::new(&BigInt::from(93750), &modulus),
                vec![x_1.pow(6), x_2.pow(25), x_3.pow(5), x_4, x_5.pow(13)],
            ),
            Term::new(
                FFE::new(&BigInt::from(112500), &modulus),
                vec![x_1.pow(6), x_2.pow(19), x_3.pow(4), x_4.pow(2), x_5.pow(21)],
            ),
            Term::new(
                FFE::new(&BigInt::from(67500), &modulus),
                vec![x_1.pow(6), x_2.pow(13), x_3.pow(3), x_4.pow(3), x_5.pow(29)],
            ),
            Term::new(
                FFE::new(&BigInt::from(20250), &modulus),
                vec![x_1.pow(6), x_2.pow(7), x_3.pow(2), x_4.pow(4), x_5.pow(37)],
            ),
            Term::new(
                FFE::new(&BigInt::from(2430), &modulus),
                vec![x_1.pow(6), x_2, x_3, x_4.pow(5), x_5.pow(45)],
            ),
            Term::new(
                FFE::new(&BigInt::from(34375), &modulus),
                vec![x_1.pow(5), x_2.pow(30), x_3.pow(5), x_4.pow(3), x_5.pow(6)],
            ),
            Term::new(
                FFE::new(&BigInt::from(103125), &modulus),
                vec![x_1.pow(5), x_2.pow(24), x_3.pow(4), x_4.pow(4), x_5.pow(14)],
            ),
            Term::new(
                FFE::new(&BigInt::from(123750), &modulus),
                vec![x_1.pow(5), x_2.pow(18), x_3.pow(3), x_4.pow(5), x_5.pow(22)],
            ),
            Term::new(
                FFE::new(&BigInt::from(74250), &modulus),
                vec![x_1.pow(5), x_2.pow(12), x_3.pow(2), x_4.pow(6), x_5.pow(30)],
            ),
            Term::new(
                FFE::new(&BigInt::from(22275), &modulus),
                vec![x_1.pow(5), x_2.pow(6), x_3, x_4.pow(7), x_5.pow(38)],
            ),
            Term::new(
                FFE::new(&BigInt::from(2673), &modulus),
                vec![x_1.pow(5), x_4.pow(8), x_5.pow(46)],
            ),
        ];
        let expected_poly = MultivariatePoly::new(terms);
        assert_eq!(poly_7.terms.len(), expected_poly.terms.len());

        assert_eq!(poly_7.terms[0], expected_poly.terms[0]);
        assert_eq!(poly_7.terms[1], expected_poly.terms[6]);
        assert_eq!(poly_7.terms[2], expected_poly.terms[1]);
        assert_eq!(poly_7.terms[3], expected_poly.terms[7]);
        assert_eq!(poly_7.terms[4], expected_poly.terms[2]);
        assert_eq!(poly_7.terms[5], expected_poly.terms[8]);
        assert_eq!(poly_7.terms[6], expected_poly.terms[3]);
        assert_eq!(poly_7.terms[7], expected_poly.terms[9]);
        assert_eq!(poly_7.terms[8], expected_poly.terms[4]);
        assert_eq!(poly_7.terms[9], expected_poly.terms[10]);
        assert_eq!(poly_7.terms[10], expected_poly.terms[5]);
        assert_eq!(poly_7.terms[11], expected_poly.terms[11]);

        // // Scenario with constant
        let seven = FFE::new(&BigInt::from(7), &modulus);
        let constant_term_seven = Term::new(seven.clone(), vec![]);
        let constant_poly = MultivariatePoly::new(vec![constant_term_seven.clone()]);

        let poly_9 = &poly_2 * &constant_poly;

        let poly_10 = &constant_poly * &poly_2;

        assert_eq!(poly_9, poly_10);

        let terms = vec![
            Term::new(
                FFE::new(&BigInt::from(175), &modulus),
                vec![x_1.pow(2), x_2.pow(12), x_3.pow(2), x_5.pow(2)],
            ),
            Term::new(
                FFE::new(&BigInt::from(210), &modulus),
                vec![x_1.pow(2), x_2.pow(6), x_3, x_4, x_5.pow(10)],
            ),
            Term::new(
                FFE::new(&BigInt::from(63), &modulus),
                vec![x_1.pow(2), x_4.pow(2), x_5.pow(18)],
            ),
        ];
        let expected_poly = MultivariatePoly::new(terms);
        assert_eq!(expected_poly, poly_9);

        let poly_11 = MultivariatePoly::new(vec![
            Term::new(FFE::one(), vec![x_1]),
            Term::new(FFE::one(), vec![x_5.pow(9)]),
            constant_term_seven,
        ]);
        let poly_12 = &poly_11 * &constant_poly;

        let poly_13 = &constant_poly * &poly_11;

        assert_eq!(poly_12, poly_13);

        let expected_poly = MultivariatePoly::new(vec![
            Term::new(FFE::new(&BigInt::from(49), &modulus), vec![]),
            Term::new(seven.clone(), vec![x_1]),
            Term::new(seven, vec![x_5.pow(9)]),
        ]);
        assert_eq!(poly_12, expected_poly)
    }

    #[test]
    fn test_add() {
        let modulus = BigInt::from(17);

        let x_1 = Var { index: 1, power: 1 };
        let x_2 = Var { index: 2, power: 1 };
        let x_3 = Var { index: 3, power: 1 };

        // SIMPLE CASE

        // x_1 + x_2 + x_3
        let term_x_1 = Term {
            coefficient: FFE::new(&BigInt::from(1), &modulus),
            vars: vec![x_1],
        };
        let poly_1 = MultivariatePoly::new(vec![term_x_1.clone()]);
        let term_x_2 = Term {
            coefficient: FFE::new(&BigInt::from(1), &modulus),
            vars: vec![x_2],
        };
        let poly_2 = MultivariatePoly::new(vec![term_x_2.clone()]);
        let term_x_3 = Term {
            coefficient: FFE::new(&BigInt::from(1), &modulus),
            vars: vec![x_3],
        };
        let poly_3 = MultivariatePoly::new(vec![term_x_3.clone()]);

        let poly_4 = &(&poly_2 + &poly_1) + &poly_3;

        let terms = vec![term_x_1, term_x_2, term_x_3];
        let expected_poly = MultivariatePoly::new(terms);
        assert_eq!(expected_poly, poly_4);

        // 3x^2 + 5y + 8z^2
        let three_x_1_square = Term {
            coefficient: FFE::new(&BigInt::from(3), &modulus),
            vars: vec![x_1.pow(2)],
        };
        let poly_5 = MultivariatePoly::new(vec![three_x_1_square.clone()]);
        let five_x_2 = Term {
            coefficient: FFE::new(&BigInt::from(5), &modulus),
            vars: vec![x_2],
        };
        let poly_6 = MultivariatePoly::new(vec![five_x_2.clone()]);
        let eight_x_3_square = Term {
            coefficient: FFE::new(&BigInt::from(8), &modulus),
            vars: vec![x_3.pow(2)],
        };
        let poly_7 = MultivariatePoly::new(vec![eight_x_3_square.clone()]);

        let poly_8 = &(&poly_5 + &poly_6) + &poly_7;

        let terms = vec![five_x_2, three_x_1_square, eight_x_3_square];
        let expected_poly = MultivariatePoly::new(terms);

        assert_eq!(expected_poly, poly_8);

        // LIKE TERMS
        let poly_9 = &poly_8 + &poly_8;

        let terms = vec![
            Term::new(FFE::new(&BigInt::from(10), &modulus), vec![x_2]),
            Term::new(FFE::new(&BigInt::from(6), &modulus), vec![x_1.pow(2)]),
            Term::new(FFE::new(&BigInt::from(16), &modulus), vec![x_3.pow(2)]),
        ];
        let expected_poly = MultivariatePoly::new(terms);
        assert_eq!(expected_poly, poly_9);
    }

    #[test]
    fn test_term_evaluate() {
        let modulus = BigInt::from(17);
        let x_1 = Var::new(1, 1);
        let x_2 = Var::new(2, 1);
        let x_3 = Var::new(3, 1);

        // x_1 + x_2^3 + x_3^5
        let term_x_1 = Term::new(FFE::one(), vec![x_1]);
        let term_x_2_exp_3 = Term::new(FFE::one(), vec![x_2.pow(3)]);
        let term_x_3_exp_5 = Term::new(FFE::one(), vec![x_3.pow(5)]);
        let poly = MultivariatePoly::new(vec![term_x_1, term_x_2_exp_3, term_x_3_exp_5]);
        let evaluation_points = vec![
            (1, FFE::new(&BigInt::from(10), &modulus)),
            (2, FFE::zero()),
            (3, FFE::new(&BigInt::from(16), &modulus)),
        ];
        let evaluation = poly.full_evaluation(&evaluation_points);
        let expected_evaluation = FFE::new(&BigInt::from(9), &modulus);
        assert_eq!(expected_evaluation, evaluation);
    }

    #[test]
    fn test_interpolate_and_full_evaluation() {
        let modulus = BigInt::from(17);
        let rounds = 1000;
        for _ in 0..rounds {
            // 2 variables
            let num_of_vars = 2;

            let num_of_points = 4;

            let mut evaluation_points: Vec<Vec<FFE>> = vec![];

            let evaluation_points = generate_random_evaluation_points(
                &mut evaluation_points,
                num_of_vars,
                num_of_points,
                &modulus,
            );

            let evaluations = generate_point(num_of_points, &modulus);

            let (poly, _, _) = MultivariatePoly::interpolate(&evaluation_points, &evaluations);

            let evaluation_point = get_evaluation_point(0, &evaluation_points);
            let evaluation = poly.full_evaluation(&evaluation_point);
            println!(
                "{:?} {:?} {:?}",
                evaluation_points, evaluations, evaluation_point
            );
            assert_eq!(evaluation, evaluations[0]);

            let evaluation_point = get_evaluation_point(1, &evaluation_points);
            let evaluation = poly.full_evaluation(&evaluation_point);
            println!(
                "{:?} {:?} {:?}",
                evaluation_points, evaluations, evaluation_point
            );
            assert_eq!(evaluation, evaluations[1]);

            let evaluation_point = get_evaluation_point(2, &evaluation_points);
            let evaluation = poly.full_evaluation(&evaluation_point);
            println!(
                "{:?} {:?} {:?}",
                evaluation_points, evaluations, evaluation_point
            );
            assert_eq!(evaluation, evaluations[2]);

            let evaluation_point = get_evaluation_point(3, &evaluation_points);
            let evaluation = poly.full_evaluation(&evaluation_point);
            assert_eq!(evaluation, evaluations[3]);
        }
        for _ in 0..rounds {
            // 3 variables
            let num_of_vars = 3;

            let num_of_points = 8;

            let mut evaluation_points: Vec<Vec<FFE>> = vec![];

            let evaluation_points = generate_random_evaluation_points(
                &mut evaluation_points,
                num_of_vars,
                num_of_points,
                &modulus,
            );

            let evaluations = generate_point(num_of_points, &modulus);

            let (poly, _, _) = MultivariatePoly::interpolate(&evaluation_points, &evaluations);

            let evaluation_point = get_evaluation_point(0, &evaluation_points);
            let evaluation = poly.full_evaluation(&evaluation_point);
            println!(
                "{:?} {:?} {:?}",
                evaluation_points, evaluations, evaluation_point
            );
            assert_eq!(evaluation, evaluations[0]);

            let evaluation_point = get_evaluation_point(1, &evaluation_points);
            let evaluation = poly.full_evaluation(&evaluation_point);
            println!(
                "{:?} {:?} {:?}",
                evaluation_points, evaluations, evaluation_point
            );
            assert_eq!(evaluation, evaluations[1]);

            let evaluation_point = get_evaluation_point(2, &evaluation_points);
            let evaluation = poly.full_evaluation(&evaluation_point);
            println!(
                "{:?} {:?} {:?}",
                evaluation_points, evaluations, evaluation_point
            );
            assert_eq!(evaluation, evaluations[2]);

            let evaluation_point = get_evaluation_point(3, &evaluation_points);
            let evaluation = poly.full_evaluation(&evaluation_point);
            println!(
                "{:?} {:?} {:?}",
                evaluation_points, evaluations, evaluation_point
            );
            assert_eq!(evaluation, evaluations[3]);

            let evaluation_point = get_evaluation_point(4, &evaluation_points);
            let evaluation = poly.full_evaluation(&evaluation_point);
            println!(
                "{:?} {:?} {:?}",
                evaluation_points, evaluations, evaluation_point
            );
            assert_eq!(evaluation, evaluations[4]);

            let evaluation_point = get_evaluation_point(5, &evaluation_points);
            let evaluation = poly.full_evaluation(&evaluation_point);
            println!(
                "{:?} {:?} {:?}",
                evaluation_points, evaluations, evaluation_point
            );
            assert_eq!(evaluation, evaluations[5]);

            let evaluation_point = get_evaluation_point(6, &evaluation_points);
            let evaluation = poly.full_evaluation(&evaluation_point);
            println!(
                "{:?} {:?} {:?}",
                evaluation_points, evaluations, evaluation_point
            );
            assert_eq!(evaluation, evaluations[6]);

            let evaluation_point = get_evaluation_point(7, &evaluation_points);
            let evaluation = poly.full_evaluation(&evaluation_point);
            println!(
                "{:?} {:?} {:?}",
                evaluation_points, evaluations, evaluation_point
            );
            assert_eq!(evaluation, evaluations[7]);
        }
    }

    #[test]
    fn test_partial_evaluation() {
        let modulus = BigInt::from(17);
        let rounds = 1000;
        for _ in 0..rounds {
            // 2 variables
            let num_of_vars = 2;

            let num_of_points = 4;

            let mut evaluation_points: Vec<Vec<FFE>> = vec![];

            let evaluation_points = generate_random_evaluation_points(
                &mut evaluation_points,
                num_of_vars,
                num_of_points,
                &modulus,
            );

            let evaluations = generate_point(num_of_points, &modulus);

            let (poly, _, _) = MultivariatePoly::interpolate(&evaluation_points, &evaluations);

            let evaluation_point = get_evaluation_point(0, &evaluation_points);
            let evaluation_1 = poly.full_evaluation(&evaluation_point);

            let first_point = vec![evaluation_point[0].clone()];
            let partial_poly_1 = poly.partial_evaluation(&first_point);
            println!("{:?} {:?}", partial_poly_1, first_point);

            let second_point = vec![evaluation_point[1].clone()];
            let evaluation_2 = partial_poly_1.full_evaluation(&second_point);

            assert_eq!(evaluation_1, evaluation_2);
        }

        for _ in 0..rounds {
            // 3 variables
            let num_of_vars = 3;

            let num_of_points = 8;

            let mut evaluation_points: Vec<Vec<FFE>> = vec![];

            let evaluation_points = generate_random_evaluation_points(
                &mut evaluation_points,
                num_of_vars,
                num_of_points,
                &modulus,
            );

            let evaluations = generate_point(num_of_points, &modulus);

            let (poly, _, _) = MultivariatePoly::interpolate(&evaluation_points, &evaluations);

            let evaluation_point = get_evaluation_point(0, &evaluation_points);
            let evaluation_1 = poly.full_evaluation(&evaluation_point);

            let first_point = vec![evaluation_point[0].clone()];
            let partial_poly_1 = poly.partial_evaluation(&first_point);
            println!("{:?} {:?}", partial_poly_1, first_point);

            let second_point = vec![evaluation_point[1].clone()];
            let partial_poly_2 = partial_poly_1.partial_evaluation(&second_point);
            println!("{:?} {:?}", partial_poly_2, second_point);

            let third_point = vec![evaluation_point[2].clone()];
            let evaluation_2 = partial_poly_2.full_evaluation(&third_point);

            assert_eq!(evaluation_1, evaluation_2);
        }
    }

    #[test]
    fn test_latex() {
        let modulus = BigInt::from(17);
        let rounds = 1000;
        for _ in 0..rounds {
            // 2 variables
            let num_of_vars = 2;

            let num_of_points = 4;

            let mut evaluation_points: Vec<Vec<FFE>> = vec![];

            let evaluation_points = generate_random_evaluation_points(
                &mut evaluation_points,
                num_of_vars,
                num_of_points,
                &modulus,
            );

            let evaluations = generate_point(num_of_points, &modulus);

            let (poly, _, _) = MultivariatePoly::interpolate(&evaluation_points, &evaluations);

            let poly_latex = poly.to_latex(num_of_vars);

            let poly_from_latex = MultivariatePoly::<FFE>::from_latex(poly_latex, &modulus);

            assert_eq!(poly, poly_from_latex);
        }

        for _ in 0..rounds {
            // 3 variables
            let num_of_vars = 3;

            let num_of_points = 8;

            let mut evaluation_points: Vec<Vec<FFE>> = vec![];

            let evaluation_points = generate_random_evaluation_points(
                &mut evaluation_points,
                num_of_vars,
                num_of_points,
                &modulus,
            );

            let evaluations = generate_point(num_of_points, &modulus);

            let (poly, _, _) = MultivariatePoly::interpolate(&evaluation_points, &evaluations);

            let poly_latex = poly.to_latex(num_of_vars);

            let poly_from_latex = MultivariatePoly::<FFE>::from_latex(poly_latex, &modulus);

            assert_eq!(poly, poly_from_latex);
        }
    }

    #[test]
    fn test_steps() {
        let modulus = BigInt::from(17);
        // let num_of_vars = 3;

        // let num_of_points = 8;

        // let mut evaluation_points: Vec<Vec<FFE>> = vec![];

        // let evaluation_points = generate_random_evaluation_points(
        //     &mut evaluation_points,
        //     num_of_vars,
        //     num_of_points,
        //     &modulus,
        // );

        // let evaluations = generate_point(num_of_points, &modulus);

        let evaluation_points = vec![
            vec![
                FFE::new(&BigInt::zero(), &modulus),
                FFE::new(&BigInt::zero(), &modulus),
                FFE::new(&BigInt::zero(), &modulus),
            ],
            vec![
                FFE::new(&BigInt::zero(), &modulus),
                FFE::new(&BigInt::zero(), &modulus),
                FFE::new(&BigInt::one(), &modulus),
            ],
            vec![
                FFE::new(&BigInt::zero(), &modulus),
                FFE::new(&BigInt::one(), &modulus),
                FFE::new(&BigInt::zero(), &modulus),
            ],
            vec![
                FFE::new(&BigInt::zero(), &modulus),
                FFE::new(&BigInt::one(), &modulus),
                FFE::new(&BigInt::one(), &modulus),
            ],
            vec![
                FFE::new(&BigInt::one(), &modulus),
                FFE::new(&BigInt::zero(), &modulus),
                FFE::new(&BigInt::zero(), &modulus),
            ],
            vec![
                FFE::new(&BigInt::one(), &modulus),
                FFE::new(&BigInt::zero(), &modulus),
                FFE::new(&BigInt::one(), &modulus),
            ],
            vec![
                FFE::new(&BigInt::one(), &modulus),
                FFE::new(&BigInt::one(), &modulus),
                FFE::new(&BigInt::zero(), &modulus),
            ],
            vec![
                FFE::new(&BigInt::one(), &modulus),
                FFE::new(&BigInt::one(), &modulus),
                FFE::new(&BigInt::one(), &modulus),
            ],
        ];
        let evaluations = vec![
            FFE::new(&BigInt::from(3), &modulus),
            FFE::new(&BigInt::from(2), &modulus),
            FFE::new(&BigInt::from(5), &modulus),
            FFE::new(&BigInt::from(7), &modulus),
            FFE::new(&BigInt::from(9), &modulus),
            FFE::new(&BigInt::zero(), &modulus),
            FFE::new(&BigInt::zero(), &modulus),
            FFE::new(&BigInt::zero(), &modulus),
        ];
        let (poly, steps, _) = MultivariatePoly::interpolate(&evaluation_points, &evaluations);
        println!("{:?}", steps);
    }
}
