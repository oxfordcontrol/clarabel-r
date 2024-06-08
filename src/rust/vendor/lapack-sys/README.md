# lapack-sys [![Package][package-img]][package-url] [![Documentation][documentation-img]][documentation-url] [![Build][build-img]][build-url]

The package provides bindings to [LAPACK] (Fortran).

## [Architecture]

## Development

The code is generated via a shell script based on the content of the `lapack`
submodule. To re-generate, run the following command:

```sh
./bin/generate.sh
```

## Contribution

Your contribution is highly appreciated. Do not hesitate to open an issue or a
pull request. Note that any contribution submitted for inclusion in the project
will be licensed according to the terms given in [LICENSE.md](LICENSE.md).

[architecture]: https://blas-lapack-rs.github.io/architecture
[lapack]: https://en.wikipedia.org/wiki/LAPACK

[build-img]: https://travis-ci.org/blas-lapack-rs/lapack-sys.svg?branch=master
[build-url]: https://travis-ci.org/blas-lapack-rs/lapack-sys
[documentation-img]: https://docs.rs/lapack-sys/badge.svg
[documentation-url]: https://docs.rs/lapack-sys
[package-img]: https://img.shields.io/crates/v/lapack-sys.svg
[package-url]: https://crates.io/crates/lapack-sys
