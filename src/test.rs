use crate::*;
use num_traits::{One, Zero};
use polynomial::Polynomial;

macro_rules! poly {
    ($($x:expr),*) => {
        Polynomial::new(vec![$(num::Rational::from_integer($x)),*])
    }
}
macro_rules! expand_poly {
    ($([$($x:expr),*]),*) => {
        vec![$(poly![$($x),*]),*].into_iter().product::<Polynomial<num::Rational>>()
    }
}

#[test]
fn test_times() {
    assert_eq!(times(42, 756), 42 * 756);
    let a = expand_poly![[2], [1, 1], [2, 1], [3, 1]];
    let c = poly![42];
    assert_eq!(times(a.clone(), 42), c * a);
}
#[test]
fn test_power() {
    type R = Polynomial<num::Rational>;
    assert_eq!(power::<i32>(3, 4), 81);
    let a = poly![1, 1];
    let b = poly![1, 4, 6, 4, 1];
    assert_eq!(power::<R>(a, 4), b);
}
#[test]
fn test_gcd() {
    assert_eq!(gcd::<i32>(0, 0), 0);
    assert_eq!(gcd::<i32>(42, 0), 42);
    assert_eq!(gcd::<i32>(0, 42), 42);
    assert_eq!(gcd::<i32>(64, 58), 2);
    assert_eq!(gcd::<i32>(97, 89), 1);
}
#[test]
fn test_gcd2() {
    type R = Polynomial<num::Rational>;
    let z = R::zero();
    assert_eq!(gcd::<R>(z.clone(), z.clone()), z.clone());
    let a = expand_poly![[2], [1, 1], [2, 1], [3, 1]];
    let b = expand_poly![[3], [1, 1], [4, 1]];
    let c = poly![1, 1];
    let d = expand_poly![[4, 1], [5, 1]];
    assert_eq!(gcd::<R>(a.clone(), z.clone()), a.clone());
    assert_eq!(gcd::<R>(z.clone(), a.clone()), a.clone());
    let mut m = gcd::<R>(a.clone(), b.clone());
    m.monic();
    assert_eq!(m, c.clone());
    let mut m = gcd::<R>(a.clone(), d.clone());
    m.monic();
    assert!(m.is_one());
}
fn check_eea<T>(a: T, b: T) -> bool
where
    T: num_traits::Zero + num_traits::One + Clone + Eq + DebugOnFeature + RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    let g = gcd::<T>(a.clone(), b.clone());
    let (d, x, y) = extended_euclidian_algorithm::<T>(a.clone(), b.clone());
    g.is_similar(&d) && &(&x * &a) + &(&y * &b) == d
}
#[test]
fn test_eea() {
    assert!(check_eea::<i32>(0, 0));
    assert!(check_eea::<i32>(42, 0));
    assert!(check_eea::<i32>(0, 42));
    assert!(check_eea::<i32>(64, 58));
    assert!(check_eea::<i32>(97, 89));
}
#[test]
fn test_eea2() {
    type R = Polynomial<num::Rational>;
    let z = R::zero();
    check_eea::<R>(z.clone(), z.clone());
    let a = expand_poly![[2], [1, 1], [2, 1], [3, 1]];
    let b = expand_poly![[3], [1, 1], [4, 1]];
    let d = expand_poly![[4, 1], [5, 1]];
    assert!(check_eea::<R>(a.clone(), z.clone()));
    assert!(check_eea::<R>(z.clone(), a.clone()));
    assert!(check_eea::<R>(a.clone(), b.clone()));
    assert!(check_eea::<R>(a.clone(), d.clone()));
}
fn check_neea<T>(a: T, b: T) -> bool
where
    T: num_traits::Zero + num_traits::One + Clone + Eq + DebugOnFeature + RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    let g = gcd::<T>(a.clone(), b.clone());
    let (d, x, y) = normalized_extended_euclidian_algorithm::<T>(a.clone(), b.clone());
    g.is_similar(&d) && &(&x * &a) + &(&y * &b) == d
}
#[test]
fn test_neea() {
    assert!(check_neea::<i32>(0, 0));
    assert!(check_neea::<i32>(42, 0));
    assert!(check_neea::<i32>(0, 42));
    assert!(check_neea::<i32>(64, 58));
    assert!(check_neea::<i32>(97, 89));
}
#[test]
fn test_neea2() {
    type R = Polynomial<num::Rational>;
    let z = R::zero();
    check_eea::<R>(z.clone(), z.clone());
    let a = expand_poly![[2], [1, 1], [2, 1], [3, 1]];
    let b = expand_poly![[3], [1, 1], [4, 1]];
    let d = expand_poly![[4, 1], [5, 1]];
    assert!(check_neea::<R>(a.clone(), z.clone()));
    assert!(check_neea::<R>(z.clone(), a.clone()));
    assert!(check_neea::<R>(a.clone(), b.clone()));
    assert!(check_neea::<R>(a.clone(), d.clone()));
}
fn check_mod_inv<T>(a: T, m: T) -> Option<T>
where
    T: num_traits::Zero + num_traits::One + Clone + Eq + DebugOnFeature + RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    modulo_inverse::<T>(a.clone(), m.clone()).map(|x| &(&(a * x) - &T::one()) % &m)
}
#[test]
fn test_mod_inv() {
    // not exists inverse
    assert_eq!(check_mod_inv::<i32>(0, 0), None);
    assert_eq!(check_mod_inv::<i32>(42, 0), None);
    assert_eq!(check_mod_inv::<i32>(0, 42), None);
    assert_eq!(check_mod_inv::<i32>(64, 58), None);
    // exists inverse
    assert_eq!(check_mod_inv::<i32>(97, 89), Some(0));
    assert_eq!(check_mod_inv::<i32>(7, 15), Some(0));
    assert_eq!(check_mod_inv::<i32>(42, 55), Some(0));
    assert_eq!(check_mod_inv::<i32>(15, 64), Some(0));
}
#[test]
fn test_mod_inv2() {
    type R = Polynomial<num::Rational>;
    // not exists inverse
    let z = R::zero();
    let a = expand_poly![[2], [1, 1], [2, 1], [3, 1]];
    let b = expand_poly![[3], [1, 1], [4, 1]];
    let d = expand_poly![[4, 1], [5, 1]];
    assert_eq!(check_mod_inv::<R>(z.clone(), z.clone()), None);
    assert_eq!(check_mod_inv::<R>(a.clone(), z.clone()), None);
    assert_eq!(check_mod_inv::<R>(z.clone(), a.clone()), None);
    assert_eq!(check_mod_inv::<R>(b.clone(), d.clone()), None);
    // exists inverse
    let sz = Some(R::zero());
    assert_eq!(check_mod_inv::<R>(a.clone(), d.clone()), sz);
    let a = poly![7, 1];
    let b = expand_poly![[3, 1], [5, 1]];
    assert_eq!(check_mod_inv::<R>(a, b), sz);
    let a = poly![42, 1];
    let b = expand_poly![[5, 1], [11, 1]];
    assert_eq!(check_mod_inv::<R>(a, b), sz);
    let a = expand_poly![[3, 1], [5, 1]];
    let b = power::<R>(poly![2, 1], 6);
    assert_eq!(check_mod_inv::<R>(a, b), sz);
}
fn check_crt<T>(u: &[T], m: &[T])
where
    T: Clone + Eq + num_traits::Zero + num_traits::One + RingNormalize + DebugOnFeature,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    let a = chinese_remainder_theorem::<T>(u, m).unwrap();
    for (u, m) in u.iter().zip(m.iter()) {
        assert!((&(&a - u) % m).is_zero());
    }
}
#[test]
fn test_crt() {
    type R = Polynomial<num::Rational>;
    let u = vec![2, 3, 2];
    let m = vec![3, 5, 7];
    check_crt::<i32>(&u, &m);
    let u = vec![
        3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4, 6, 2, 6, 4, 3, 3,
    ];
    let m = vec![
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97, 101,
    ];
    check_crt::<i128>(&u, &m);
    let u = vec![
        poly![3, 1],
        poly![5],
        poly![7, 1],
        poly![1, 1],
        poly![2],
        poly![1, 1],
        poly![3, 3],
        poly![1, 7],
    ];
    let m = vec![
        poly![1, 1, 1],
        poly![1, 2, 1],
        poly![1, 3, 1],
        poly![1, 4, 1],
        poly![1, 5, 1],
        poly![1, 6, 1],
        poly![1, 7, 1],
        poly![1, 8, 1],
    ];
    check_crt::<R>(&u, &m);
}
