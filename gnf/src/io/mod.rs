extern crate quick_error;

pub mod refflat;

pub mod error {

    use std::io::Error as StdIoError;
    use csv::Error as CsvError;
    use FeatureError;

    quick_error! {
        #[derive(Debug)]
        pub enum ParseError {
            Column(err: CsvError) {
                description(err.description())
                from()
            }
            Io(err: StdIoError) {
                description(err.description())
                from()
            }
            FormatSpecific(string: &'static str) {
                description(string)
            }
            Feature(err: FeatureError) {
                description(err.description())
                from()
            }
        }
    }
}
