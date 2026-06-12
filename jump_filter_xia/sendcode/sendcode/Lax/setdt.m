function dt = setdt(msh, U, t, tc, quad, bs, cfl, tp)

% Evaluate means of density, momentum, energy, pressure and sound speed
% at each element
um = 0.5 * reshape(quad.weights' * (bs.phi{1} * U), [3, msh.nLElems]);
pm = computePressure(um(1, :), um(2, :), um(3, :), tc);
cm = computeSpeedOfSound(um(1, :), pm, tc);

% dissipation coefficient
alpha = max(abs(um(2, :) ./ um(1, :)) + cm);

% Set the time step
dt = cfl * min(msh.elemLength(:, msh.LElems)) / ((alpha - 1) * (alpha > 1.0e-9) + 1);
if (t + dt) > tp
    dt = tp - t;
end

end

