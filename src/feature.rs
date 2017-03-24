//! Interval-based annotation features.

use bio::utils::{Interval, IntervalError};


pub trait NamedInterval: Sized {

    /// Underlying interval struct.
    fn interval(&self) -> &Interval<u64>;

    /// Name of the interval.
    fn name(&self) -> Option<&str>;

    /// Name setter that returns the implementor itself.
    ///
    /// This function is expected to mutate the implementing type.
    fn with_name<T: Into<String>>(self, name: T) -> Self;

    /// Coordinate setter that returns the implementor itself.
    ///
    /// This function is expected to mutate the implementing type.
    fn with_coords(self, start: u64, end: u64) -> Result<Self, IntervalError>;

    /// Start coordinate of the interval.
    fn start(&self) -> u64 {
        self.interval().start
    }

    /// End coordinate of the interval.
    fn end(&self) -> u64 {
        self.interval().end
    }

    /// The number of bases covered by the interval.
    fn span(&self) -> u64 {
        self.end() - self.start()
    }

    /// Whether two intervals have an intersection or not.
    fn overlaps(&self, other: &Self) -> bool {
        (other.start() <= self.start() && self.start() <= other.end()) ||
        (other.start() <= self.end() && self.end() <= other.end())
    }

    /// Whether one interval completely contains the other.
    fn envelops(&self, other: &Self) -> bool {
        self.start() <= other.start() && self.end() >= other.end()
    }

    /// Whether two intervals cover a contiguous region without any overlaps.
    fn adjacent(&self, other: &Self) -> bool {
        self.end() == other.start() || self.start() == other.end()
    }
}

/// Macro for default function implementations of interval types.
macro_rules! impl_ninterval {
    ($struct_ty:ty) => (

        impl NamedInterval for $struct_ty {

            /// Name of the interval.
            fn name(&self) -> Option<&str> {
                self.name.as_ref().map(|n| n.as_str())
            }

            fn with_name<T>(mut self, name: T) -> $struct_ty
                where T: Into<String>
            {
                self.name = Some(name.into());
                self
            }

            fn interval(&self) -> &Interval<u64> {
                &self.interval
            }

            fn with_coords(mut self, start: u64, end: u64) -> Result<$struct_ty, IntervalError> {
                Interval::new(start..end)
                    .map(|iv| {
                        self.interval = iv;
                        self
                    })
            }
        }

    );
}

/// Default implementation of the `Interval` trait.
///
/// This struct also provides static methods for creating exons, transcripts, and genes.
#[derive(Debug)]
pub struct Feature {
    interval: Interval<u64>,
    name: Option<String>,
}

impl Default for Feature {

    fn default() -> Feature {
        Feature { interval: Interval::new(0..0).unwrap(), name: None }
    }
}

impl_ninterval!(Feature);

#[cfg(test)]
mod test_feature {
    use super::*;

    #[test]
    fn default() {
        let fx = Feature::default();
        assert_eq!(fx.start(), 0);
        assert_eq!(fx.end(), 0);
        assert_eq!(fx.name(), None);
    }

    #[test]
    fn with_name() {
        let fx1 = Feature::default()
            .with_name("fx1");
        assert_eq!(fx1.start(), 0);
        assert_eq!(fx1.end(), 0);
        assert_eq!(fx1.name(), Some("fx1"));

        let fx2 = Feature::default()
            .with_name("fx2".to_owned());
        assert_eq!(fx2.start(), 0);
        assert_eq!(fx2.end(), 0);
        assert_eq!(fx2.name(), Some("fx2"));
    }

    #[test]
    fn with_coords() {
        let fxm = Feature::default()
            .with_coords(1, 3);
        assert!(fxm.is_ok());
        let fx = fxm.unwrap();
        assert_eq!(fx.start(), 1);
        assert_eq!(fx.end(), 3);
        assert_eq!(fx.name(), None);
    }

    #[test]
    fn with_coords_err() {
        let fxm = Feature::default()
            .with_coords(3, 1);
        assert!(fxm.is_err());
    }

    #[test]
    fn with_multiples() {
        let fxm = Feature::default()
            .with_coords(20, 30)
            .map(|f| f.with_name("fx"));
        assert!(fxm.is_ok());
        let fx = fxm.unwrap();
        assert_eq!(fx.start(), 20);
        assert_eq!(fx.end(), 30);
        assert_eq!(fx.name(), Some("fx"));
    }
}
