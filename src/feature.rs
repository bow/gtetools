//! Interval-based annotation features.

use bio::utils::{Interval, IntervalError};


pub trait NamedInterval: Sized {

    /// Underlying interval struct.
    fn interval(&self) -> &Interval<u64>;  // TODO: Generalize over interval types.

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

    /// Whether two intervals have an overlap or not.
    fn overlaps(&self, other: &Self) -> bool {
        self.start() < other.end() && other.start() < self.end()
    }

    /// Whether one interval completely contains the other.
    fn envelops(&self, other: &Self) -> bool {
        self.start() <= other.start() && self.end() >= other.end()
    }

    /// Whether two intervals cover a contiguous region without any overlaps.
    fn is_adjacent(&self, other: &Self) -> bool {
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

    fn make_feature(start: u64, end: u64) -> Feature {
        Feature::default().with_coords(start, end).unwrap()
    }

    #[test]
    fn span() {
        let fx = make_feature(0, 15);
        assert_eq!(fx.span(), 15);
    }

    #[test]
    fn overlaps() {
        let fx1 = make_feature(100, 115);

        let fx2 = make_feature(110, 120);
        assert!(fx1.overlaps(&fx2));
        assert!(fx1.overlaps(&fx2));

        let fx3 = make_feature(115, 120);
        assert!(!fx1.overlaps(&fx3));
        assert!(!fx3.overlaps(&fx1));

        let fx4 = make_feature(90, 100);
        assert!(!fx1.overlaps(&fx4));
        assert!(!fx4.overlaps(&fx1));

        let fx5 = make_feature(200, 300);
        assert!(!fx1.overlaps(&fx5));
        assert!(!fx5.overlaps(&fx1));
    }

    #[test]
    fn envelops() {
        let fx1 = make_feature(100, 120);

        let fx2 = make_feature(105, 115);
        assert!(fx1.envelops(&fx2));
        assert!(!fx2.envelops(&fx1));

        let fx3 = make_feature(100, 105);
        assert!(fx1.envelops(&fx3));
        assert!(!fx3.envelops(&fx1));

        let fx4 = make_feature(115, 120);
        assert!(fx1.envelops(&fx4));
        assert!(!fx4.envelops(&fx1));

        let fx5 = make_feature(90, 105);
        assert!(!fx1.envelops(&fx5));
        assert!(!fx5.envelops(&fx1));

        let fx6 = make_feature(115, 130);
        assert!(!fx1.envelops(&fx5));
        assert!(!fx6.envelops(&fx1));

        let fx7 = make_feature(80, 90);
        assert!(!fx1.envelops(&fx7));
        assert!(!fx7.envelops(&fx1));
    }

    #[test]
    fn is_adjacent() {
        let fx1 = make_feature(100, 120);

        let fx2 = make_feature(90, 100);
        assert!(fx1.is_adjacent(&fx2));
        assert!(fx2.is_adjacent(&fx1));

        let fx3 = make_feature(120, 130);
        assert!(fx1.is_adjacent(&fx3));
        assert!(fx3.is_adjacent(&fx1));

        let fx4 = make_feature(90, 99);
        assert!(!fx1.is_adjacent(&fx4));
        assert!(!fx4.is_adjacent(&fx1));

        let fx5 = make_feature(119, 130);
        assert!(!fx1.is_adjacent(&fx5));
        assert!(!fx5.is_adjacent(&fx1));

        let fx6 = make_feature(100, 110);
        assert!(!fx1.is_adjacent(&fx6));
        assert!(!fx6.is_adjacent(&fx1));

        let fx7 = make_feature(110, 120);
        assert!(!fx1.is_adjacent(&fx7));
        assert!(!fx7.is_adjacent(&fx1));
    }
}
