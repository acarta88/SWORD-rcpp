# Contributing to SWORD

Thank you for your interest in contributing to the SWORD package.

## Reporting bugs

Please open an issue at <https://github.com/acarta88/SWORD-rcpp/issues> with:

- A minimal reproducible example
- Your R version and OS
- The output of `sessionInfo()`

## Suggesting enhancements

Open an issue describing the proposed feature and the use case that motivates it.

## Pull requests

1. Fork the repository and create a branch from `main`.
2. Make your changes and ensure `R CMD check` passes with no errors or warnings.
3. Add or update tests in `tests/testthat/` as appropriate.
4. Submit the pull request with a clear description of the changes.

## Code style

- Follow the existing code conventions in `R/SWORD.R`.
- Document all exported functions with `roxygen2`.
- C++ code lives in `src/`; follow the existing Rcpp conventions.

## Contact

Andrea Carta — <acarta88@yahoo.it>
