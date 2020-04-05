#![cfg_attr(feature = "__internal_inject_debug", recursion_limit = "8")]
mod sealed {
    pub trait SizedExt: std::marker::Sized + std::fmt::Debug + std::fmt::Display {}
    impl<T> SizedExt for T where T: std::marker::Sized + std::fmt::Debug + std::fmt::Display {}
    #[cfg(not(feature = "__internal_inject_debug"))]
    pub use std::marker::Sized;
    #[cfg(feature = "__internal_inject_debug")]
    pub use SizedExt as Sized;
}
mod ring_traits;
#[cfg(test)]
mod test;
pub use ring_traits::{EuclideanRingOperation, RingNormalize};

/** calcurate $`pa`$ with mutliprecation by doubling
```
use ring_algorithm::times;
assert_eq!(times::<i32>(2, 16), 32);
```
*/
pub fn times<T>(a: T, mut p: u64) -> T
where
    T: sealed::Sized + num_traits::Zero + for<'x> std::ops::AddAssign<&'x T>,
    for<'x> &'x T: std::ops::Add<Output = T>,
{
    let mut x = T::zero();
    let mut y = a;
    loop {
        if p % 2 == 1 {
            x += &y;
        }
        p /= 2;
        if p == 0 {
            break;
        }
        y = &y + &y;
    }
    x
}

/** calcurate $`a^p`$ with exponentiation by squaring
```
use ring_algorithm::power;
assert_eq!(power::<i32>(2, 16), 65536);
```
*/
pub fn power<T>(a: T, mut p: u64) -> T
where
    T: sealed::Sized + num_traits::One + for<'x> std::ops::MulAssign<&'x T>,
    for<'x> &'x T: std::ops::Mul<Output = T>,
{
    let mut x = T::one();
    let mut y = a;
    loop {
        if p % 2 == 1 {
            x *= &y;
        }
        p /= 2;
        if p == 0 {
            break;
        }
        y = &y * &y;
    }
    x
}

/** calcurate greatest common divisor
```
use ring_algorithm::gcd;
assert_eq!(gcd::<i32>(15, 21), 3);
assert_eq!(gcd::<i32>(14, 15), 1);
assert_eq!(gcd::<i32>(0, 42), 42);
assert_eq!(gcd::<i32>(0, 0), 0);
```
*/
pub fn gcd<T>(mut x: T, mut y: T) -> T
where
    T: sealed::Sized + num_traits::Zero,
    for<'x> &'x T: std::ops::Rem<Output = T>,
{
    while !y.is_zero() {
        let r = &x % &y;
        x = y;
        y = r;
    }
    x
}

/** test $`\gcd(x, y) = 1`$
*/
pub fn is_coprime<T>(x: T, y: T) -> bool
where
    T: sealed::Sized + Eq + num_traits::Zero + num_traits::One + RingNormalize,
    for<'x> &'x T: std::ops::Rem<Output = T>,
{
    gcd::<T>(x, y).into_normalize().is_one()
}

/** extended euclidian algorithm

calcurate g (`gcd(a, b)`) and x, y ( $`g = ax + by`$ )
```
use ring_algorithm::{gcd, extended_euclidian_algorithm};
let a = 314;
let b = 271;
let (d, x, y) = extended_euclidian_algorithm::<i32>(a, b);
assert_eq!(d, gcd::<i32>(a, b));
assert_eq!(d, x * a + y * b);
```
 */
pub fn extended_euclidian_algorithm<T>(x: T, y: T) -> (T, T, T)
where
    T: sealed::Sized + num_traits::Zero + num_traits::One,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    let mut old = (x, T::one(), T::zero());
    let mut now = (y, T::zero(), T::one());
    while !now.0.is_zero() {
        let q = &old.0 / &now.0;
        let new = (
            &old.0 - &(&q * &now.0),
            &old.1 - &(&q * &now.1),
            &old.2 - &(&q * &now.2),
        );
        old = now;
        now = new;
    }
    old
}

/** extended euclidian algorithm with normalize
```
use ring_algorithm::{gcd, normalized_extended_euclidian_algorithm, RingNormalize};
let a = 314;
let b = 271;
let (d, x, y) = normalized_extended_euclidian_algorithm::<i32>(a, b);
assert_eq!(d, gcd::<i32>(a, b));
assert_eq!(d, x * a + y * b);
```
*/
pub fn normalized_extended_euclidian_algorithm<T>(x: T, y: T) -> (T, T, T)
where
    T: sealed::Sized + num_traits::Zero + num_traits::One + RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    let lc_x = x.leading_unit();
    let lc_y = y.leading_unit();
    let mut old = (x.into_normalize(), &T::one() / &lc_x, T::zero());
    let mut now = (y.into_normalize(), T::zero(), &T::one() / &lc_y);
    while !now.0.is_zero() {
        let q = &old.0 / &now.0;
        let r = &old.0 % &now.0;
        let lc_r = r.leading_unit();
        let new = (
            r.into_normalize(),
            &(&old.1 - &(&q * &now.1)) / &lc_r,
            &(&old.2 - &(&q * &now.2)) / &lc_r,
        );
        old = now;
        now = new;
    }
    old
}

/** calc inverse in modulo

calc x ($`ax \equiv 1 \pmod{m}`$)
```
use ring_algorithm::modulo_inverse;
let a = 42;
let m = 55;
let b = modulo_inverse::<i32>(a, m).unwrap();
assert_eq!((a * b - 1) % m, 0);
```
*/
pub fn modulo_inverse<T>(a: T, m: T) -> Option<T>
where
    T: sealed::Sized + Eq + num_traits::Zero + num_traits::One + RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    let (gcd, inv_a, _) = normalized_extended_euclidian_algorithm::<T>(a, m);
    if gcd.is_one() {
        Some(inv_a)
    } else {
        None
    }
}

/** division in modulo

calc x ($`bx \equiv a \pmod{m}`$)
```
use ring_algorithm::modulo_divison;
let a = 42;
let b = 32;
let m = 98;
let x = modulo_divison::<i32>(a, b, m).unwrap();
assert_eq!((b * x - a) % m, 0);
```
*/
pub fn modulo_divison<T>(a: T, b: T, m: T) -> Option<T>
where
    T: sealed::Sized + Eq + num_traits::Zero + num_traits::One + RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    let (gcd, inv_b, _) = normalized_extended_euclidian_algorithm::<T>(b, m);
    if (&a % &gcd).is_zero() {
        Some(&a / &gcd * inv_b)
    } else {
        None
    }
}

/** Chinese remainder theorem

```
use ring_algorithm::chinese_remainder_theorem;
let u = vec![2, 3, 2];
let m = vec![3, 5, 7];
let a = chinese_remainder_theorem::<i32>(&u, &m).unwrap();
for (u, m) in u.iter().zip(m.iter()) {
    assert_eq!((a - u) % m, 0);
}
```
*/
pub fn chinese_remainder_theorem<T>(u: &[T], m: &[T]) -> Option<T>
where
    T: sealed::Sized + Clone + Eq + num_traits::Zero + num_traits::One + RingNormalize,
    for<'x> &'x T: EuclideanRingOperation<T>,
{
    if u.len() != m.len() {
        return None;
    }
    let mut v = Vec::with_capacity(u.len());
    for (i, (u_i, m_i)) in u.iter().zip(m.iter()).enumerate() {
        let coef_i = modulo_inverse::<T>(
            m[0..i].iter().fold(T::one(), |p, v| &(&p * v) % m_i),
            m_i.clone(),
        )?;
        let t = v
            .iter()
            .zip(m.iter())
            .rev()
            .fold(T::zero(), |t, (v_j, m_j)| &(&(m_j * &t) + v_j) % m_i);
        v.push(&(&(u_i - &t) * &coef_i) % m_i);
    }
    let mut ret = v.pop().unwrap();
    for (v_i, m_i) in v.iter().zip(m.iter()).rev() {
        ret = &(&ret * m_i) + v_i;
    }
    Some(ret)
}
