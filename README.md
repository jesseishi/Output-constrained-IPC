# Output-constrained invididual pitch control

This README will describe the code found in this repository entry and describe how this code
and data was used towards the publication "Output-constrained individual pitch control
methods using the multiblade coordinate transformation: Trading off actuation effort and
blade fatigue load reduction for wind turbines". [Input the proper citation once a
preprint is available.]

Data and code for WES publication: "Output-constrained individual pitch control
methods using the multiblade coordinate transformation: Trading off actuation effort and
blade fatigue load reduction for wind turbines" [![DOI](https://data.4tu.nl/v3/datasets/372325a3-306e-4578-9c72-4fcda690a999/doi-badge.svg)](https://doi.org/10.4121/372325a3-306e-4578-9c72-4fcda690a999)

In case you have any questions, please contact me:
[j.i.s.hummel@tudelft.nl](mailto:j.i.s.hummel@tudelft.nl).

## Getting started

You can get started in two ways:

1. Define a case in `generate_cases.py`, run the case with `Run_cIPC.m` or multiple
   cases with `Run_batch_cIPC.m`, and then analyze with a script inspired by
   `Plot_WES.m` or `Plot_NAWEA.m`.
2. Download our results and put them in the `Results/data` folder. You can now analyze
   all the results using `Plot_WES.m`.

## Folder structure

- `IEA-15-240-RWT`
  - Defines the IEA 15 MW reference turbine. Adjusted from
    [GitHub](https://github.com/IEAWindTask37/IEA-15-240-RWT), see
    `IEA-15-240-RWT/README.md` for the changes made.
- `preplot-postplot` (optional)
  - My plotting utilities to make nice plots with less boilerplate, available on [GitHub](https://github.com/jesseishi/preplot-postplot).
- `Results` (optional)
  - `data`: Stores simulation data to be used by plotting scripts. If you want to
      recreate our results, store our datasets in this folder. This folder has
      folders that mirror the folders of `InputFiles`. Then for each folder,
      there are folders indicating what reference signal was used by the controllers
      (e.g. `ref0` for a 0 Nm reference), then finally for each simulation there is an
      `.outb` file with OpenFAST outputs and a `.csv` files which also contrains logged
      signals from Simulink.
  - `figures`: Stores resulting figures.
  - `InputFiles`: Stores input files used by OpenFAST to run the different cases.
- `src`:
  - `Helper`: Contains helper functions.
  - `Lib`: Library with some static data (linearizations and optimal azimuth offsets).
  - `Models`: Constains the different Simulink models.
  - `TempCache` (optimal): Temporary cache folder.
  - `generate_cases.py`: Generate input files to run different cases (wind speed,
turbulence intensity, etc...).
  - `Plot_XYZ.m`: Scripts to plot the results.
  - `README.md`: This README.
  - `Run_batch_cIPC.m`: Script to run many OpenFAST simulations with ℓ2 and ℓ∞
      output-constrained IPC.
  - `TurbSim_x64.exe`: TurbSim executable to generate turbulent wind files when running
    `generate_cases.py`.

## Software used

- Python 3.12.1 with the [OpenFAST toolbox](https://github.com/OpenFAST/python-toolbox)
- Matlab and Simulink R2022a with the
  [preplot-postplot](https://github.com/jesseishi/preplot-postplot) toolbox for plotting
- OpenFAST 3.5.0
- TurbSim v2.0

## Licence

Apache-2.0 license, see `LICENSE`.
