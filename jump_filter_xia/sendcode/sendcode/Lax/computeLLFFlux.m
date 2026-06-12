% Compute local Lax-Friedrichs flux
function H = computeLLFFlux(UL, UR, pl, pr, FL, FR, tc)

% normal velocity (inward or outward)
ul = UL(2, :) ./ UL(1, :);
ur = UR(2, :) ./ UR(1, :);

% speed of sound
cl = computeSpeedOfSound(UL(1, :), pl, tc);
cr = computeSpeedOfSound(UR(1, :), pr, tc);

% dissipation coefficient
alpha = max(abs(ul) + cl, abs(ur) + cr);

% local Lax-Friedrichs flux
H = 0.5 * (FL + FR - alpha .* (UR - UL));

end