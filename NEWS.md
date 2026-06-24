# Rpcop 1.2.2

- Added rOpenSci-oriented package review infrastructure, including pkgcheck,
  R-CMD-check, coverage, and pkgdown workflows.
- Added SRR standards documentation, package citation metadata, a pkgdown
  configuration, contributor guidance, and a PCOP vignette.
- Hardened `pcop()` input validation for unsupported types, non-finite data,
  tuning-parameter ranges, and plotting preconditions.
- Returned `pcop` S3 objects while preserving the original `pcop.f1` and
  `pcop.f2` list elements, with new `print()`, `summary()`, and `plot()`
  methods.
