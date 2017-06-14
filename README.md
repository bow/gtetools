# gtetools

[![Build Status](https://travis-ci.org/bow/gtetools.svg?branch=master)](https://travis-ci.org/bow/gtetools)

`gtetools` is a suite of command line tools for working with various gene annotation formats.


# Installation

Since gtetools is still in alpha, please follow the [development setup](#development-setup) instructions below.


# Development Setup

gtetools is written in the Rust programming language (v1.18 or above). If you do not have the Rust toolchain
installed yet, you will need to install it. The Rust [official installation guide](https://www.rustup.rs/) is
easy to follow and should get you going in no time.

If you already have toolchain installed, ensure that you are at least on version 1.18 and execute the following
commands:

    # Clone the repository and cd into it
    $ git clone https://github.com/bow/gtetools.git
    $ cd gtetools

    # Install the dependencies
    $ cargo update

    # Run the test suite (if you want to play around with the code)
    $ cargo test

    # Build the release version (to ./target/release/gtetools)
    $ cargo build --release


# License

gtetools is licensed under the BSD license.
