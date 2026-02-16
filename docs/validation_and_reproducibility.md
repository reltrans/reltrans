# Validation and Reproducibility

## Numerical Validation

The repository includes numerical validation benchmarks in the `Benchmarks/` directory.

These benchmarks compare model outputs against reference data to ensure
numerical consistency across builds and environments.

They are intended to verify scientific correctness rather than measure runtime performance.

---

## Running Validation Benchmarks

Please refer to the instructions provided within the `Benchmarks/`
directory for building and executing the validation program.

---

## Reproducibility Considerations

For consistent scientific results, users may consider:

- Keeping grid resolution fixed when comparing outputs
- Recording all physical parameter values used
- Using consistent compiler settings
- Documenting floating-point precision
- Noting hardware details when comparing performance

---

## Validation vs Performance

The existing benchmarks focus on validating numerical results.

Performance profiling and runtime scaling analysis are separate concerns
and may vary depending on system configuration.

