% the flux function 
function F = computeF(rho, m, E, p)

F = [m; m.^2 ./ rho + p; m ./ rho .* (E + p)];

end

