# AutoIR: Automated Infrared Spectra from Molecular Dynamics

AutoIR is a Python package for computing infrared (IR) spectra from molecular dynamics (MD) simulations. Starting from a SMILES string, it automates the entire workflow:

1. **System setup**: generates molecular topologies and simulation boxes for gas phase, pure liquids, and multi-component mixtures using the OpenFF toolkit.
2. **MD simulation**: runs NVT/NPT molecular dynamics with OpenMM, recording dipole moments along the trajectory.
3. **Spectral analysis**: computes the IR spectrum via the dipole autocorrelation function and its Fourier transform, including correction factors.

## Installation

### From repository

Clone the repository and install with `pip`:

```bash
git clone https://github.com/digital-synthesis-lab/autoir.git
cd autoir
pip install .
```

## Usage

### Command line

Once installed, the `autoir` command provides subcommands for different simulation types.

#### Gas phase

Run a gas-phase simulation for a single molecule given as a SMILES string:

```bash
autoir gas "CCO"
```

Common options:

| Option | Default | Description |
|--------|---------|-------------|
| `-o`, `--output` | auto UUID | Output directory |
| `--temperature` | 300 | Temperature (K) |
| `--time_step` | 2 | Timestep (fs) |
| `--nvt_equi_steps` | 100,000 | Equilibration steps |
| `--nvt_prod_steps` | 300,000 | Production steps |
| `--trj_freq` | 10,000 | Trajectory dump frequency (steps) |
| `--platform` | auto | OpenMM platform (e.g. `CUDA`, `CPU`) |
| `--seed` | 12345 | Random seed |

#### Liquid phase

Run an NPT + NVT simulation of a pure liquid:

```bash
autoir liquid "CCO" -n 80 --temperature 300 --pressure 1.0
```

Additional options beyond the gas phase:

| Option | Default | Description |
|--------|---------|-------------|
| `-n`, `--n_mols` | 80 | Number of molecules in the box |
| `--pressure` | 1.0 | Pressure (atm) |
| `--target_density` | 0.5 | Initial packing density (g/cm³) |
| `--npt_equi_steps` | 500,000 | NPT equilibration steps |
| `--nvt_equi_steps` | 100,000 | NVT equilibration steps |
| `--nvt_prod_steps` | 500,000 | NVT production steps |

#### Binary mixture

Simulate a liquid mixture of two components given their SMILES strings and a molar ratio:

```bash
autoir mixture "CCO" "O" --ratio 0.3 -n 100
```

| Option | Default | Description |
|--------|---------|-------------|
| `--ratio` | 0.5 | Molar fraction of the first component |

All other options from the `liquid` command are also available.

#### N-component mixture

For mixtures with more than two components, use `n_mixture`:

```bash
autoir n_mixture -s "CCO" -s "O" -s "CC" -r "0.5 0.3 0.2"
```

| Option | Default | Description |
|--------|---------|-------------|
| `-s`, `--smiles` | — | SMILES string (repeat for each component) |
| `-r`, `--ratio` | uniform | Molar ratios (repeat or space-separated; normalized automatically) |

All other liquid-phase options apply.

#### Analyzing simulation output

After a simulation, compute the IR spectrum from the recorded dipole moment file:

```bash
autoir analyze -i dipoles.csv -o ir.csv
```

| Option | Default | Description |
|--------|---------|-------------|
| `-i`, `--input` | `dipoles.csv` | Dipole moment time series from the MD run |
| `-o`, `--output` | `ir.csv` | Output IR spectrum file |
| `-t`, `--truncate_autocorr` | 5000 | Number of steps used for the autocorrelation function |

#### Single-point energy

Compute the potential energy for a single frame of a PDB file (for example, a
structure produced by a prior simulation). One SMILES string must be supplied
for each unique molecule present in the PDB, since OpenFF cannot infer
bond orders from PDB coordinates alone:

```bash
autoir singlepoint trajectory.pdb -s "CCO" --frame 0
```

| Option | Default | Description |
|--------|---------|-------------|
| `-s`, `--smiles` | — | SMILES for each unique molecule in the PDB (repeatable) |
| `--frame` | 0 | Zero-based MODEL index to evaluate |
| `--force_field` | `openff_unconstrained-2.0.0.offxml` | SMIRNOFF force-field filename |
| `--platform` | auto | OpenMM platform (e.g. `CUDA`, `CPU`) |

The command prints the potential energy in kJ/mol to stdout. The periodic box
is taken from the PDB `CRYST1` record when present; otherwise a 4 nm cubic
non-periodic box is used (matching the gas-phase workflow).

#### Getting help

```bash
autoir --help
autoir gas --help
autoir liquid --help
autoir mixture --help
autoir n_mixture --help
autoir analyze --help
autoir singlepoint --help
```

### Output files

Each simulation creates a directory (specified with `-o`, or a random UUID by default) containing:

| File | Description |
|------|-------------|
| `dipoles.csv` | Dipole moment (x, y, z) recorded every `trj_freq` steps |
| `equi.csv` | Thermodynamic data during equilibration |
| `prod.csv` | Thermodynamic data during production |
| `avgE.csv` | Running-average potential energy |
| `trajectory.pdb` | Atomic trajectory |
| `ir.csv` | Computed IR spectrum (wavenumber vs. intensity) |
| `job.json` | Simulation parameters and summary metrics |

### API

#### Gas phase simulation

```python
from autoir.create import gas_simulation
from autoir.openmm.simulator import OpenMMSimulator
from autoir.analyze import process_file, smooth_ir

topology, interchange = gas_simulation("CCO")

sim = OpenMMSimulator(temperature=300, pressure=0.0, nvt_prod_steps=300_000)
sim.run(interchange, npt_equi_steps=0, nvt_equi_steps=100_000, nvt_prod_steps=300_000)

df = process_file("dipoles.csv", out_file="ir.csv")
wn, ir = smooth_ir(df)
```

#### Liquid phase simulation

```python
from autoir.create import liq_simulation
from autoir.openmm.simulator import OpenMMSimulator
from autoir.analyze import process_file, smooth_ir

topology, interchange = liq_simulation("CCO", n_mols=80, target_density=0.5)

sim = OpenMMSimulator(temperature=300, pressure=1.0)
sim.run(interchange, npt_equi_steps=500_000, nvt_equi_steps=100_000, nvt_prod_steps=500_000)

df = process_file("dipoles.csv", out_file="ir.csv")
wn, ir = smooth_ir(df)
```

#### Binary mixture simulation

```python
from autoir.create import mix_simulation
from autoir.openmm.simulator import OpenMMSimulator
from autoir.analyze import process_file

topology, interchange = mix_simulation("CCO", "O", ratio=0.3, n_mols=100)

sim = OpenMMSimulator(temperature=300, pressure=1.0)
sim.run(interchange, npt_equi_steps=500_000, nvt_equi_steps=100_000, nvt_prod_steps=500_000)

df = process_file("dipoles.csv", out_file="ir.csv")
```

#### N-component mixture simulation

```python
from autoir.create import n_mix_simulation
from autoir.openmm.simulator import OpenMMSimulator
from autoir.analyze import process_file

smiles = ["CCO", "O", "CC"]
ratios = [0.5, 0.3, 0.2]

topology, interchange = n_mix_simulation(smiles, ratio_list=ratios, n_mols=100)

sim = OpenMMSimulator(temperature=300, pressure=1.0)
sim.run(interchange, npt_equi_steps=500_000, nvt_equi_steps=100_000, nvt_prod_steps=500_000)

df = process_file("dipoles.csv", out_file="ir.csv")
```

#### Computing the IR spectrum from a dipole file

```python
from autoir.analyze import compute_autocorr, compute_spectra, smooth_ir
import pandas as pd

autocorr, timestep = compute_autocorr("dipoles.csv", truncate_autocorr=5000)
wavenums, spectra = compute_spectra(autocorr, timestep, T=300)

wn, ir = smooth_ir(pd.DataFrame({"w": wavenums, "IR": spectra}))
```

#### Single-point energy from a PDB frame

```python
from autoir.singlepoint import single_point_energy

energy_kj_per_mol = single_point_energy(
    pdb_file="trajectory.pdb",
    smiles=["CCO"],   # one SMILES per unique molecule in the PDB
    frame=0,
)
```

## Partial charges and OpenEye

AutoIR parameterizes molecules with the OpenFF Toolkit, which chooses the partial
charge method based on what is installed in your environment. If the OpenEye
Toolkits are installed and a valid license is found, OpenFF assigns AM1BCC-ELF10
charges, averaged over a set of selected conformers. If OpenEye is not available,
OpenFF falls back to single-conformer AM1BCC charges from AmberTools with a
warning.

To use the OpenEye Toolkits set the `OE_LICENSE`
environment variable to the path of your license file. Academic licenses are
available from [OpenEye](https://www.eyesopen.com/academic-licensing). 


If you use OpenEye with AutoIR, please also cite

OpenEye Toolkits <VERSION>. OpenEye, Cadence Molecular Sciences, Santa Fe, NM.
http://www.eyesopen.com.


## Citing

If you use AutoIR in a publication, please cite the following paper:

```bibtex
@article{melle2026irid,
  author = {Melle, Yannah J. U. and Nguyen, Thanh and Lopez, Jeffrey and Schwalbe-Koda, Daniel},
  title = {Automatic identification of compounds in molecular mixtures from liquid-phase infrared spectra},
  journal = {Chem. Sci.},
  year = {2026},
  publisher = {The Royal Society of Chemistry},
  doi = {10.1039/D6SC01583B},
  url = {http://dx.doi.org/10.1039/D6SC01583B},
}
```

## License

AutoIR is distributed under the BSD-3-Clause license.

SPDX: BSD-3-Clause



