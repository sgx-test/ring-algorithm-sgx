use crate::ring_traits::{DebugOnFeature, RingNormalize};
use num_traits::{One, Zero};
use std::ops::*;

#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Polynomial<T> {
    coef: Vec<T>,
}

macro_rules! from_assign {
    ($t:ident, $op_t:ident, $op:ident, $assign:ident, $require:tt) => {
        impl<$t: $require<$t>> $op_t for Polynomial<$t> {
            type Output = Self;
            fn $op(mut self, other: Self) -> Self::Output {
                self.$assign(&other);
                self
            }
        }
        impl<'a, $t: $require<$t>> $op_t for &'a Polynomial<$t> {
            type Output = Polynomial<$t>;
            fn $op(self, other: Self) -> Self::Output {
                let mut t = self.clone();
                t.$assign(other);
                t
            }
        }
    };
}

// additive monoid
pub trait AddAssignRequire<M>: Clone + Zero + for<'x> AddAssign<&'x M> + DebugOnFeature {}
impl<M> AddAssignRequire<M> for M where M: Clone + Zero + for<'x> AddAssign<&'x M> + DebugOnFeature {}
impl<'a, M: AddAssignRequire<M>> AddAssign<&'a Polynomial<M>> for Polynomial<M> {
    fn add_assign(&mut self, other: &Polynomial<M>) {
        let len = self.len();
        self.extend(other.len());
        self.coef
            .iter_mut()
            .zip(other.coef.iter())
            .for_each(|(l, r)| *l += r);
        if len == other.len() {
            self.trim_zero()
        }
    }
}
from_assign!(M, Add, add, add_assign, AddAssignRequire);
impl<M: AddAssignRequire<M>> Zero for Polynomial<M> {
    fn zero() -> Self {
        Self { coef: Vec::new() }
    }
    fn is_zero(&self) -> bool {
        self.deg().is_none()
    }
}
impl<M: AddAssignRequire<M>> std::iter::Sum for Polynomial<M> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::zero(), Add::add)
    }
}
impl<M: DebugOnFeature> Polynomial<M> {
    fn len(&self) -> usize {
        self.coef.len()
    }
    pub fn deg(&self) -> Option<usize> {
        if self.coef.is_empty() {
            None
        } else {
            Some(self.len() - 1)
        }
    }
    pub fn lc(&self) -> Option<&M> {
        self.deg().map(|d| &self.coef[d])
    }
    fn trim_zero(&mut self)
    where
        M: Zero,
    {
        let len = self
            .coef
            .iter()
            .rposition(|x| !x.is_zero())
            .map(|pos| pos + 1)
            .unwrap_or(0);
        self.coef.truncate(len);
    }
    pub fn new(coef: Vec<M>) -> Self
    where
        M: Zero,
    {
        let mut poly = Self { coef };
        poly.trim_zero();
        poly
    }
    fn extend(&mut self, len: usize)
    where
        M: Clone + Zero,
    {
        if len > self.len() {
            self.coef.resize(len, M::zero());
        }
    }
    /// c*x^d (c=coefficent, d=degree)
    pub fn from_monomial(coefficent: M, degree: usize) -> Self
    where
        M: Clone + Zero,
    {
        let mut coef = vec![M::zero(); degree + 1];
        coef[degree] = coefficent;
        Self { coef }
    }
}

// additive group
impl<G> Neg for Polynomial<G>
where
    G: Neg<Output = G> + DebugOnFeature,
{
    type Output = Self;
    fn neg(self) -> Self {
        Polynomial {
            coef: self.coef.into_iter().map(|v| -v).collect(),
        }
    }
}
pub trait SubAssignRequire<G>: Clone + Zero + for<'x> SubAssign<&'x G> + DebugOnFeature {}
impl<G> SubAssignRequire<G> for G where G: Clone + Zero + for<'x> SubAssign<&'x G> + DebugOnFeature {}
impl<'a, G: SubAssignRequire<G>> SubAssign<&'a Polynomial<G>> for Polynomial<G> {
    fn sub_assign(&mut self, other: &Polynomial<G>) {
        let len = self.len();
        self.extend(other.len());
        self.coef
            .iter_mut()
            .zip(other.coef.iter())
            .for_each(|(l, r)| *l -= r);
        if len == other.len() {
            self.trim_zero()
        }
    }
}
from_assign!(G, Sub, sub, sub_assign, SubAssignRequire);

// unitary ring
fn mul_aux<R>(sum: &mut [R], coef: &R, vec: &[R])
where
    R: AddAssignRequire<R>,
    for<'x> &'x R: Mul<Output = R>,
{
    sum.iter_mut()
        .zip(vec.iter())
        .for_each(|(l, r)| *l += &(coef * r));
}
impl<'a, R> Mul for &'a Polynomial<R>
where
    R: AddAssignRequire<R>,
    for<'x> &'x R: Mul<Output = R>,
{
    type Output = Polynomial<R>;
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn mul(self, other: Self) -> Self::Output {
        if self.is_zero() || other.is_zero() {
            return Self::Output::zero();
        }
        let mut coef = vec![R::zero(); self.len() + other.len() - 1];
        self.coef
            .iter()
            .enumerate()
            .for_each(|(i, c)| mul_aux::<R>(&mut coef[i..], c, &other.coef));
        Polynomial::<R>::new(coef) // R may not be a domain.
    }
}
impl<R> Mul for Polynomial<R>
where
    R: AddAssignRequire<R>,
    for<'x> &'x R: Mul<Output = R>,
{
    type Output = Self;
    fn mul(self, other: Self) -> Self::Output {
        &self * &other
    }
}
impl<'a, R> MulAssign<&'a Polynomial<R>> for Polynomial<R>
where
    R: AddAssignRequire<R>,
    for<'x> &'x R: Mul<Output = R>,
{
    fn mul_assign(&mut self, other: &Polynomial<R>) {
        *self = &*self * other;
    }
}
impl<R> One for Polynomial<R>
where
    R: AddAssignRequire<R> + One,
    for<'x> &'x R: Mul<Output = R>,
{
    fn one() -> Self {
        Self {
            coef: vec![R::one()],
        }
    }
}
impl<R> std::iter::Product for Polynomial<R>
where
    R: AddAssignRequire<R> + One,
    for<'x> &'x R: Mul<Output = R>,
{
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::one(), Mul::mul)
    }
}
impl<R> std::fmt::Display for Polynomial<R>
where
    R: std::cmp::Eq + std::fmt::Display + Zero + One,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let vec = &self.coef;
        if vec.is_empty() {
            return write!(f, "0");
        }
        let mut is_first = true;
        for i in (0..vec.len()).rev() {
            let c = &vec[i];
            if c.is_zero() {
                continue;
            }
            if is_first {
                is_first = false;
            } else {
                write!(f, "+")?
            }
            if c.is_one() {
                match i {
                    0 => write!(f, "1")?,
                    1 => write!(f, "x")?,
                    _ => write!(f, "x^{}", i)?,
                }
            } else {
                match i {
                    0 => write!(f, "{}", c)?,
                    1 => write!(f, "{}*x", c)?,
                    _ => write!(f, "{}*x^{}", c, i)?,
                }
            }
        }
        Ok(())
    }
}
impl<R: DebugOnFeature> Polynomial<R> {
    pub fn eval<'a>(&self, x: &'a R) -> R
    where
        R: AddAssignRequire<R> + One + MulAssign<&'a R>,
        for<'x> &'x R: Mul<Output = R>,
    {
        if self.coef.is_empty() {
            return R::zero();
        }
        let mut sum = self.coef[self.len() - 1].clone();
        for i in (0..self.len() - 1).rev() {
            sum *= x;
            sum += &self.coef[i];
        }
        sum
    }
}

// division
impl<'a, K> Div for &'a Polynomial<K>
where
    K: AddAssignRequire<K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    type Output = Polynomial<K>;
    fn div(self, other: Self) -> Self::Output {
        let mut f = self.clone();
        f.division(other)
    }
}
impl<'a, K> DivAssign<&'a Polynomial<K>> for Polynomial<K>
where
    K: AddAssignRequire<K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    fn div_assign(&mut self, other: &Polynomial<K>) {
        *self = &*self / other;
    }
}
impl<'a, K> RemAssign<&'a Polynomial<K>> for Polynomial<K>
where
    K: AddAssignRequire<K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    fn rem_assign(&mut self, other: &Polynomial<K>) {
        self.division(other);
    }
}
impl<'a, K> Rem for &'a Polynomial<K>
where
    K: AddAssignRequire<K> + for<'x> SubAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
{
    type Output = Polynomial<K>;
    fn rem(self, other: Self) -> Self::Output {
        let mut t = self.clone();
        t %= other;
        t
    }
}
impl<K> RingNormalize for Polynomial<K>
where
    K: AddAssignRequire<K> + One + for<'x> DivAssign<&'x K>,
    for<'x> &'x K: Mul<Output = K>,
{
    fn leading_unit(&self) -> Self {
        if let Some(x) = self.lc() {
            Self::from_monomial(x.clone(), 0)
        } else {
            Self::one()
        }
    }
    fn normalize_mut(&mut self) {
        self.monic();
    }
}
impl<K: DebugOnFeature> Polynomial<K> {
    pub fn monic(&mut self)
    where
        K: Clone + for<'x> DivAssign<&'x K>,
    {
        if let Some(lc) = self.lc() {
            let lc = lc.clone();
            self.coef.iter_mut().for_each(|v| *v /= &lc);
        }
    }
    pub fn division(&mut self, other: &Self) -> Self
    where
        K: AddAssignRequire<K> + for<'x> SubAssign<&'x K>,
        for<'x> &'x K: Mul<Output = K> + Div<Output = K>,
    {
        let g_deg = other.deg().expect("Division by zero");
        if self.deg() < other.deg() {
            return Self::zero();
        }
        let lc = other.lc().unwrap();
        let mut coef = vec![K::zero(); self.len() - other.len() + 1];
        while self.deg() >= other.deg() {
            let d = self.deg().unwrap() - g_deg;
            let c = self.lc().unwrap() / lc;
            for i in 0..other.len() - 1 {
                self.coef[i + d] -= &(&c * &other.coef[i]);
            }
            self.coef.pop(); // new deg < prev deg
            self.trim_zero();
            coef[d] = c;
        }
        Self { coef }
    }
}
