% Monte Carlo Simulation of the 1D Ising Model - Quick Version
% This is a faster version with reduced parameters for quick demonstration

clear all;
close all;
clc;

%% Parameters (reduced for faster execution)
N = 50;                     % Number of spins in the chain
J = 1.0;                    % Coupling constant (positive for ferromagnetic)
T_min = 0.1;                % Minimum temperature
T_max = 5.0;                % Maximum temperature
T_steps = 20;               % Number of temperature points
equilibration_steps = 2000; % Steps to reach equilibrium
measurement_steps = 3000;   % Steps for measurements

%% Temperature range
T = linspace(T_min, T_max, T_steps);

%% Storage for results
magnetization = zeros(1, T_steps);
magnetization_std = zeros(1, T_steps);
energy = zeros(1, T_steps);

%% Main simulation loop
fprintf('Starting 1D Ising Model Simulation (Quick Version)\n');
fprintf('System size: %d spins\n', N);
fprintf('Temperature range: %.2f to %.2f\n\n', T_min, T_max);

for t_idx = 1:T_steps
    T_current = T(t_idx);
    beta = 1.0 / T_current;  % Inverse temperature

    % Initialize spins randomly
    spins = 2 * randi([0, 1], 1, N) - 1;  % Random +1 or -1

    % Equilibration phase
    for step = 1:equilibration_steps
        spins = metropolis_step(spins, N, J, beta);
    end

    % Measurement phase
    mag_samples = zeros(1, measurement_steps);
    energy_samples = zeros(1, measurement_steps);

    for step = 1:measurement_steps
        spins = metropolis_step(spins, N, J, beta);

        % Calculate magnetization per spin
        mag_samples(step) = abs(sum(spins)) / N;

        % Calculate energy per spin
        energy_samples(step) = calculate_energy(spins, N, J) / N;
    end

    % Store average values
    magnetization(t_idx) = mean(mag_samples);
    magnetization_std(t_idx) = std(mag_samples);
    energy(t_idx) = mean(energy_samples);

    % Progress update
    if mod(t_idx, 5) == 0
        fprintf('Progress: %d/%d temperatures completed\n', t_idx, T_steps);
    end
end

fprintf('\nSimulation complete!\n\n');

%% Plotting
figure('visible', 'off');
set(gcf, 'Position', [100, 100, 1200, 500]);

% Plot magnetization vs temperature
subplot(1, 2, 1);
errorbar(T, magnetization, magnetization_std, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Temperature (k_BT/J)', 'FontSize', 12);
ylabel('Magnetization |M|', 'FontSize', 12);
title('1D Ising Model: Magnetization vs Temperature', 'FontSize', 14);
grid on;
ylim([0, 1.1]);
set(gca, 'FontSize', 11);

% Plot energy vs temperature
subplot(1, 2, 2);
plot(T, energy, 'r-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Temperature (k_BT/J)', 'FontSize', 12);
ylabel('Energy per Spin (E/N)', 'FontSize', 12);
title('1D Ising Model: Energy vs Temperature', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 11);

% Save the figure
print('ising_model_results_quick.png', '-dpng', '-r150');
fprintf('Results saved to ising_model_results_quick.png\n');
