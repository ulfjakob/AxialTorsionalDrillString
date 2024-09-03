# AxialTorsionalDrillString
Code for the numerical experiments in "Axial and torsional self-excited vibrations of a distributed drill-string" JSV 2019
https://www.sciencedirect.com/science/article/pii/S0022460X1830854X

This repository contains MATLAB code for the analysis and simulation of axial-torsional vibrations in drill-string systems. The code models the stability and dynamic behavior of drill-strings under various conditions, focusing on the bit-rock interaction and boundary conditions.

## Directory Structure

Each folder contains the code to run the corresponding experiments that are shown in the paper.

- **Axial**: Model with only the axial dynamics in the drill string, assuming constand torsional velocity
- **Axial - nondim**: Only axial dynamics in the non-dimensional formulation
- **Coupled**: Model with both axial and torsional dynamics
- **Coupled - nondim**: Axial-torsional model in the non-dim formulation

## Key Files

- `MAIN.m`: Each experiment is run by executing a MAIN_[experimentName] file
- `parameters.m`: Simulation and parameters are specified here..
- `stabilityMap.m`: Generates stability maps.

## License

This project is licensed under the MIT License.

## Usage

Run simulations by executing `MAIN_[experimentName].m` in MATLAB within the desired directory.
