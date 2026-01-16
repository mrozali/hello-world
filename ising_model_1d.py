#!/usr/bin/env python3
"""
Monte Carlo Simulation of the 1D Ising Model
This script simulates the 1D Ising model using the Metropolis algorithm
and plots the magnetization as a function of temperature
"""

import numpy as np
import matplotlib.pyplot as plt

def metropolis_step(spins, N, J, beta):
    """Perform one Monte Carlo sweep using Metropolis algorithm"""
    for _ in range(N):
        # Select a random spin
        site = np.random.randint(0, N)

        # Calculate energy change if this spin is flipped
        # For 1D with periodic boundary conditions
        left = (site - 1) % N
        right = (site + 1) % N

        # Energy change = 2 * J * spin[site] * (spin[left] + spin[right])
        dE = 2 * J * spins[site] * (spins[left] + spins[right])

        # Metropolis acceptance criterion
        if dE <= 0 or np.random.random() < np.exp(-beta * dE):
            spins[site] *= -1

    return spins

def calculate_energy(spins, N, J):
    """Calculate total energy of the system"""
    E = 0
    for i in range(N):
        right = (i + 1) % N
        E -= J * spins[i] * spins[right]
    return E

# Parameters
N = 50                      # Number of spins in the chain
J = 1.0                     # Coupling constant (positive for ferromagnetic)
T_min = 0.1                 # Minimum temperature
T_max = 5.0                 # Maximum temperature
T_steps = 25                # Number of temperature points
equilibration_steps = 2000  # Steps to reach equilibrium
measurement_steps = 3000    # Steps for measurements

# Temperature range
T = np.linspace(T_min, T_max, T_steps)

# Storage for results
magnetization = np.zeros(T_steps)
magnetization_std = np.zeros(T_steps)
energy = np.zeros(T_steps)

print('Starting 1D Ising Model Simulation')
print(f'System size: {N} spins')
print(f'Temperature range: {T_min:.2f} to {T_max:.2f}\n')

# Main simulation loop
for t_idx, T_current in enumerate(T):
    beta = 1.0 / T_current  # Inverse temperature

    # Initialize spins randomly
    spins = np.random.choice([-1, 1], size=N)

    # Equilibration phase
    for _ in range(equilibration_steps):
        spins = metropolis_step(spins, N, J, beta)

    # Measurement phase
    mag_samples = np.zeros(measurement_steps)
    energy_samples = np.zeros(measurement_steps)

    for step in range(measurement_steps):
        spins = metropolis_step(spins, N, J, beta)

        # Calculate magnetization per spin
        mag_samples[step] = np.abs(np.sum(spins)) / N

        # Calculate energy per spin
        energy_samples[step] = calculate_energy(spins, N, J) / N

    # Store average values
    magnetization[t_idx] = np.mean(mag_samples)
    magnetization_std[t_idx] = np.std(mag_samples)
    energy[t_idx] = np.mean(energy_samples)

    # Progress update
    if (t_idx + 1) % 5 == 0:
        print(f'Progress: {t_idx + 1}/{T_steps} temperatures completed')

print('\nSimulation complete!\n')

# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot magnetization vs temperature
ax1.errorbar(T, magnetization, yerr=magnetization_std,
             fmt='b-o', linewidth=2, markersize=6, capsize=3)
ax1.set_xlabel('Temperature (k$_B$T/J)', fontsize=12)
ax1.set_ylabel('Magnetization |M|', fontsize=12)
ax1.set_title('1D Ising Model: Magnetization vs Temperature', fontsize=14)
ax1.grid(True, alpha=0.3)
ax1.set_ylim([0, 1.1])
ax1.tick_params(labelsize=11)

# Plot energy vs temperature
ax2.plot(T, energy, 'r-o', linewidth=2, markersize=6)
ax2.set_xlabel('Temperature (k$_B$T/J)', fontsize=12)
ax2.set_ylabel('Energy per Spin (E/N)', fontsize=12)
ax2.set_title('1D Ising Model: Energy vs Temperature', fontsize=14)
ax2.grid(True, alpha=0.3)
ax2.tick_params(labelsize=11)

plt.tight_layout()
plt.savefig('ising_model_results.png', dpi=150, bbox_inches='tight')
print('Results saved to ising_model_results.png')
plt.close()
