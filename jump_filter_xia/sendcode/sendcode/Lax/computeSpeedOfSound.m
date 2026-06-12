% speed of sound
function c = computeSpeedOfSound(rho, p, tc)

c = sqrt(abs(tc.gamma * p ./ rho));

end