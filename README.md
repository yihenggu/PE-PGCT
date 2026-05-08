# PE-PGCT

This repository provides R code for the power-enhanced panel Granger causality test used in the paper.

## Files

- `granger_utils_PEGranger.R`: main utility functions, including `PE_PGCT()`.
- `example_data.csv`: example panel dataset for demonstration.
- `main_example.R`: example showing how to run the test.

## Required R packages

Please install the following packages before running the example:

- `plm`
- `lmtest`
- `harmonicmeanp`

## Notes

- The example dataset is simulated and is intended only to demonstrate usage.
- The input data for `PE_PGCT()` should be a balanced panel with columns `id`, `time`, `x`, and `y`.
