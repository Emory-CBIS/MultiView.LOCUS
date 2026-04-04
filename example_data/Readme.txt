Example data overview
=====================

This folder provides the datasets used in the examples.

Main data files
---------------

These are the files you will usually use directly:

- **Simulation I**
  - `Data_x_simu1.rds`
  - `Data_y_simu1.rds`

- **Simulation II**
  - `Data_x_simu2.rds`
  - `Data_y_simu2.rds`

- **Simulation III**
  - `Data_x.rds`
  - `Data_y.rds`

The info below is  about how we get the Simulation data above,  optional
---------------------------------------------------------------------------------------------------

The `S/` folder stores source objects used to generate the simulation data.

- `S_xreal.rds`, `S_yreal.rds`
- `simulated_common_SourceMar7.RData`
- `simulatedspecificSourceMar6.RData`
- `simulated_sources_type2.Rdata`

How `data_generating.R` is used
-------------------------------

- `data_generating.R` reads source objects from `S/`.
- It generates the `Data_x*` and `Data_y*` output files listed above.
