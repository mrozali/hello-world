function E = calculate_energy(spins, N, J)
    % Calculate total energy of the system
    % E = -J * sum over nearest neighbors of s_i * s_j

    E = 0;
    for i = 1:N
        right = mod(i, N) + 1;  % Right neighbor with periodic BC
        E = E - J * spins(i) * spins(right);
    end
end
