# README - Final Project for Heterogeneous Agent Model Course

This repository contains the code for the final project in the Heterogeneous Agent Model course. It includes scripts to generate Impulse Response Functions (IRFs) for an MIT shock under various scenarios.

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

Considering that the majority of the code is sample code provided in class, you might want to keep the repository private to avoid potential issues with academic integrity and copyright. Always check with your course instructor or institution's policies regarding sharing course materials online.
