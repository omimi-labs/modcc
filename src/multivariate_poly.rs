use crate::ff::FiniteFieldElement;

#[derive(Debug, Default, Clone)]
struct Term {
    index: usize,
    power: usize,
}

#[derive(Debug, Default, Clone)]
pub struct CombinedTerms<F> {
    coefficient: F,
    terms: Vec<Term>,
}

#[derive(Debug, Clone)]
pub struct MultivariatePoly<F> {
    terms: Vec<CombinedTerms<F>>,
}

pub trait Polynomial<F>: Sized {
    fn is_zero(&self) -> bool;

    fn is_one(&self) -> bool;

    fn new(terms: CombinedTerms<F>) -> Self;

    fn zero() -> Self;

    fn one() -> Self;

    fn terms(&self) -> &Vec<CombinedTerms<F>>;
}

impl<F: FiniteFieldElement> Polynomial<F> for MultivariatePoly<F> {
    fn is_zero(&self) -> bool {}

    fn is_one(&self) -> bool {}

    fn new(terms: CombinedTerms<F>) -> Self {}

    fn one() -> Self {}

    fn zero() -> Self {
        let zero_combined_term = CombinedTerms::default();
        MultivariatePoly {
            terms: vec![zero_combined_term],
        }
    }

    fn terms(&self) -> &Vec<CombinedTerms<F>> {
        todo!()
    }
}
