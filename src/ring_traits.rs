pub trait UnitaryRing<Output = Self>:
    Sized
    + Clone
    + num_traits::Zero
    + num_traits::One
    + std::ops::Neg<Output = Output>
    + for<'x> std::ops::AddAssign<&'x Self>
    + for<'x> std::ops::SubAssign<&'x Self>
    + for<'x> std::ops::MulAssign<&'x Self>
{
}
impl<T> UnitaryRing for T where
    T: Sized
        + Clone
        + num_traits::Zero
        + num_traits::One
        + std::ops::Neg<Output = Self>
        + for<'x> std::ops::AddAssign<&'x Self>
        + for<'x> std::ops::SubAssign<&'x Self>
        + for<'x> std::ops::MulAssign<&'x Self>
{
}
pub trait EuclideanRing<Output = Self>:
    UnitaryRing + for<'x> std::ops::DivAssign<&'x Self> + for<'x> std::ops::RemAssign<&'x Self>
{
}
impl<T> EuclideanRing for T where
    T: UnitaryRing + for<'x> std::ops::DivAssign<&'x Self> + for<'x> std::ops::RemAssign<&'x Self>
{
}
pub trait RingOperation<Output = Self>:
    Sized
    + std::ops::Add<Output = Output>
    + std::ops::Sub<Output = Output>
    + std::ops::Mul<Output = Output>
{
}
impl<T> RingOperation<T> for T where
    T: Sized + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + std::ops::Mul<Output = T>
{
}
impl<'a, T> RingOperation<T> for &'a T where
    &'a T:
        Sized + std::ops::Add<Output = T> + std::ops::Sub<Output = T> + std::ops::Mul<Output = T>
{
}
pub trait EuclideanRingOperation<Output = Self>:
    RingOperation<Output> + std::ops::Div<Output = Output> + std::ops::Rem<Output = Output>
{
}
impl<T> EuclideanRingOperation<T> for T where
    T: RingOperation<T> + std::ops::Div<Output = T> + std::ops::Rem<Output = T>
{
}
impl<'a, T> EuclideanRingOperation<T> for &'a T where
    &'a T: RingOperation<T> + std::ops::Div<Output = T> + std::ops::Rem<Output = T>
{
}
/** Normarize ring element

`abs(a)` in $`\mathbb{Z}`$.
`a/lc(a)` in $`R[x]`$ (`lc(x)` is leading coefficent of x).
*/
pub trait RingNormalize {
    fn leading_unit(&self) -> Self;
    fn normalize_mut(&mut self);
    fn into_normalize(mut self) -> Self
    where
        Self: Sized,
    {
        self.normalize_mut();
        self
    }
    fn normalize(&self) -> Self
    where
        Self: Clone,
    {
        self.clone().into_normalize()
    }
    fn is_similar(&self, other: &Self) -> bool
    where
        Self: Clone + Eq,
    {
        self.normalize() == other.normalize()
    }
}

macro_rules! ring_normalize {
    ($t:ty) => {
        impl RingNormalize for $t {
            fn leading_unit(&self) -> Self {
                use num_traits::{One, Zero};
                if self >= &Self::zero() {
                    Self::one()
                } else {
                    -Self::one()
                }
            }
            fn normalize_mut(&mut self) {
                *self = self.abs();
            }
        }
    };
}

ring_normalize!(i8);
ring_normalize!(i16);
ring_normalize!(i32);
ring_normalize!(i64);
ring_normalize!(i128);
ring_normalize!(isize);

impl RingNormalize for num_bigint::BigInt {
    fn leading_unit(&self) -> Self {
        use num_traits::One;
        if self.sign() == num_bigint::Sign::Minus {
            -Self::one()
        } else {
            Self::one()
        }
    }
    fn normalize_mut(&mut self) {
        if self.sign() == num_bigint::Sign::Minus {
            *self = -&*self
        }
    }
}
