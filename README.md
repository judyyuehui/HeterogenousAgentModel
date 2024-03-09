# README - Computational Project for ECON 33540

This repository contains the code for my computational project for Macroeconomics with Heterogeneous Households. It includes scripts to generate IRFs for various types of MIT shocks.

## Getting Started

To obtain the IRFs for different types of shocks, follow these instructions:

- **TFP Shock:**
  - Run `ge_ayagari_irf.m`
  - Set `param.tfpshock_size = -0.05`
  - Keep other `shocksize` values at 0
  - Change `plot_path` to your local directory

- **Discount Rate Shock:**
  - Run `ge_ayagari_irf.m`
  - Set `param.discountshock_size = -0.001`
  - Keep other `shocksize` values at 0
  - Change `plot_path` to your local directory

- **Preference Shock:**
  - Run `ge_ayagari_irf.m`
  - Set `param.preferenceshock_size = 0.001`
  - Keep other `shocksize` values at 0
  - Change `plot_path` to your local directory

- **Borrowing Constraint (approximated with asset quality):**
  - Run `ge_ayagari_irf_borrowlimit2.m`
  - Change `plot_path` to your local directory

## Notes

- Ensure that all parameter changes are reverted back to their original values after running each script to avoid conflicts.

## Repository Privacy

Considering that the majority of the code is sample code provided in class, please let me know if I should keep the repository private to avoid potential issues with academic integrity and copyright. 
