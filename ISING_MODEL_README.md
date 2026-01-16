# 1D Ising Model Monte Carlo Simulation

This MATLAB script simulates the one-dimensional Ising model using the Metropolis Monte Carlo algorithm.

## Overview

The Ising model is a mathematical model of ferromagnetism in statistical mechanics. In the 1D case, spins are arranged in a chain where each spin interacts only with its nearest neighbors.

### Hamiltonian

The energy of the system is given by:

```
H = -J Σ s_i * s_{i+1}
```

where:
- `J` is the coupling constant (positive for ferromagnetic interaction)
- `s_i` is the spin at site i (either +1 or -1)
- The sum is over all nearest-neighbor pairs

## Features

- **Metropolis Algorithm**: Efficient Monte Carlo sampling using the Metropolis acceptance criterion
- **Periodic Boundary Conditions**: The chain forms a ring (first and last spins are neighbors)
- **Temperature Sweep**: Automatically simulates over a range of temperatures
- **Statistical Analysis**: Computes average magnetization with error bars
- **Visualization**: Generates plots of magnetization and energy vs temperature

## Usage

Simply run the script in MATLAB:

```matlab
ising_model_1d
```

## Parameters

You can modify these parameters at the top of the script:

- `N`: Number of spins (default: 100)
- `J`: Coupling constant (default: 1.0)
- `T_min`, `T_max`: Temperature range (default: 0.1 to 5.0)
- `T_steps`: Number of temperature points to sample (default: 30)
- `equilibration_steps`: Steps to reach thermal equilibrium (default: 5000)
- `measurement_steps`: Steps for taking measurements (default: 10000)

## Output

The script produces:

1. **Console output**: Progress updates during simulation
2. **Figure window**: Two plots showing:
   - Magnetization vs Temperature (with error bars)
   - Energy per spin vs Temperature
3. **Saved file**: `ising_model_results.png` containing the plots

## Physical Insights

### 1D Ising Model Behavior

Unlike the 2D Ising model, the 1D Ising model:

- **No phase transition** at finite temperature (only at T = 0)
- Magnetization decreases smoothly with temperature
- At T = 0: Perfect alignment (M = 1)
- As T → ∞: Random orientation (M → 0)

This is consistent with the **Peierls argument** and the **Mermin-Wagner theorem**, which state that continuous symmetries cannot be spontaneously broken at finite temperature in systems with sufficiently short-range interactions in dimensions d ≤ 2.

## Algorithm Details

### Metropolis Algorithm

For each Monte Carlo step:

1. Select a random spin
2. Calculate the energy change ΔE if the spin is flipped
3. Accept the flip if:
   - ΔE ≤ 0 (energy decreases), or
   - With probability exp(-βΔE) (Boltzmann factor)
4. Repeat N times per sweep

### Measurements

- **Equilibration phase**: System reaches thermal equilibrium
- **Measurement phase**: Statistical averages are computed
- **Magnetization**: Average of |Σ s_i| / N over all measurements

## Requirements

- MATLAB (any recent version)
- No additional toolboxes required

## Typical Run Time

- Approximately 10-30 seconds depending on your system

## Extension Ideas

1. Vary system size N to study finite-size effects
2. Compute specific heat and magnetic susceptibility
3. Implement different boundary conditions (open vs periodic)
4. Add external magnetic field
5. Study autocorrelation times and critical slowing down
