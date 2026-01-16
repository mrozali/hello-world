function spins = metropolis_step(spins, N, J, beta)
    % Perform one Monte Carlo sweep using Metropolis algorithm
    % Each sweep attempts N spin flips

    for i = 1:N
        % Select a random spin
        site = randi(N);

        % Calculate energy change if this spin is flipped
        % For 1D with periodic boundary conditions
        left = mod(site - 2, N) + 1;   % Left neighbor
        right = mod(site, N) + 1;       % Right neighbor

        % Energy change = 2 * J * spin[i] * (spin[left] + spin[right])
        dE = 2 * J * spins(site) * (spins(left) + spins(right));

        % Metropolis acceptance criterion
        if dE <= 0
            % Always accept if energy decreases
            spins(site) = -spins(site);
        else
            % Accept with probability exp(-beta * dE)
            if rand() < exp(-beta * dE)
                spins(site) = -spins(site);
            end
        end
    end
end
