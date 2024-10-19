# Motion After Effect (MAE)

## Overview

This repository contains two MATLAB simulations related to neural dynamics: `MAE_Simulations` and `BinocularRivalryTanh`. These simulations model the behavior of neural populations using differential equations and nonlinear functions. These programs were used to simulate the Motion After Effect (MAE) model presented in the paper [Sustained Rhythmic Brain Activity Underlies Visual Motion Perception in Zebrafish](https://www.cell.com/cell-reports/fulltext/S2211-1247(16)31321-3).

## Programs

1. **MAE_Simulations**
   - Simulates a model with direction-selective populations and comparator populations.
   - Implements a nonlinear function to saturate firing rates.
   - Uses a Runge-Kutta method for integration.

2. **BinocularRivalryTanh**
   - Models binocular rivalry using a tanh activation function.
   - Simulates two populations with adaptation variables and incorporates noise.
   - Uses a 6-step Runge-Kutta Fehlberg method for integration.

## Requirements

- MATLAB (preferably version 2016 or later).
- No additional toolboxes are required.

## Usage

### MAE_Simulations

To run the `MAE_Simulations` program, simply call the function from the MATLAB command window:

```matlab
MAE_Simulations
```

This will generate and save simulation data in a file named `File_to_save.mat`.

### BinocularRivalryTanh

To run the `BinocularRivalryTanh` program, use the following command:

```matlab
BinocularRivalryTanh
```

This will produce plots for firing rates and adaptation variables for the two populations.

## Functions

- **MAE_Simulations**
  - Defines global variables for time and input currents.
  - Uses a differential equation solver (`HH`) to update state variables.
  - Includes a saturation function (`f`) to constrain firing rates.

- **BinocularRivalryTanh**
  - Similar structure with a tanh activation function.
  - Incorporates noise and pulse excitations to the populations.
  - Also uses a differential equation solver (`HH`) to update state variables.

## Output

- For `MAE_Simulations`: A `.mat` file containing the variables `u1` and `u2`.
- For `BinocularRivalryTanh`: Plots of firing rates and adaptation variables.

## License

This code is provided for educational purposes. Please credit the authors if used in publications or presentations.

