# r-src : A Source of `BLAS` and `LAPACK` via R

[![Github Actions Build Status](https://github.com/blas-lapack-rs/r-src/workflows/Tests/badge.svg)](https://github.com/blas-lapack-rs/r-src/actions)
[![Crates.io](https://img.shields.io/crates/v/r-src.svg)](https://crates.io/crates/r_src)
[![Documentation](https://docs.rs/extendr-api/badge.svg)](https://docs.rs/r_src)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Note that only double precision routines are provided by `R`. The build
script for this particular crate owes a lot to
[libR-sys](https://github.com/extendr/libR-sys) of the
[extendr](https://github.com/extendr) project, with sections lifted
from there, in addition to new code.

Use of this package is explained
[here](https://blas-lapack-rs.github.io/usage), specifically under the
_Sources_ section.

## Installation

This library can be installed as usual. It can also built from source
but has extra platform-specific dependencies for `Windows` such as
`msys2`, typically provided by
[`Rtools`](https://cran.r-project.org/bin/windows/Rtools/).

## Configuration

Two components are required to build the library:

1. [`R`](https://cran.r-project.org/): It needs to be installed and
   available in the search path. On `Windows`, the `Rtools` binaries
   are also expected to be on the search path.
2. [`Rust`](https://www.rust-lang.org/learn/get-started): It is
   recommended to install `Rust` using `rustup`; search path should
   include `Rust` binaries.

Once `R` and `Rust` are configured, the library can be easily built:

```bash
# macOS & Linux
cargo build

# Windows
cargo build --target x86_64-pc-windows-gnu
```

To test the build, run `cargo test`.

```bash
# macOS & Linux
cargo test

# Windows
cargo test --target x86_64-pc-windows-gnu
```

## Contribution

Your contribution is highly appreciated. Do not hesitate to open an issue or a
pull request. Note that any contribution submitted for inclusion in the project
will be licensed according to the terms given in [LICENSE.md](LICENSE.md).

[architecture]: https://blas-lapack-rs.github.io/architecture
[blas]: https://en.wikipedia.org/wiki/BLAS
[lapack]: https://en.wikipedia.org/wiki/LAPACK

[build-img]: https://travis-ci.org/blas-lapack-rs/r-src.svg?branch=master
[build-url]: https://travis-ci.org/blas-lapack-rs/r-src
[documentation-img]: https://docs.rs/r_src/badge.svg
[documentation-url]: https://docs.rs/r_src
[package-img]: https://img.shields.io/crates/v/r_src.svg
[package-url]: https://crates.io/crates/r_src


