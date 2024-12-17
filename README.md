# autoir

The `autoir` package is designed to facilitate the analysis and simulation of infrared spectra using molecular dynamics simulations. It provides tools for setting up simulations, rendering input files, and analyzing the results to compute infrared spectra.

## Modules

- **analyze.py**: Contains functions to compute autocorrelation functions and spectra from simulation data.
- **const.py**: Defines physical constants used throughout the package.
- **create.py**: Provides functions to set up molecular simulations, including estimating box sizes and generating topologies for gas and liquid phases.
- **render.py**: Handles the rendering of LAMMPS input files using Jinja2 templates.
- **lammps/template.in**: A template file for LAMMPS simulations, which is populated with job-specific parameters.

## CLI

The `autoir.cli` module contains scripts to implement a command line interface for the `autoir` package, allowing users to interact with the package's functionality from the terminal.
