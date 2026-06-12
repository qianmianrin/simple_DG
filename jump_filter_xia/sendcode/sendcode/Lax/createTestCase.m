function tc = createTestCase(gamma, pro)

% ratio of specific heats, with gamma = 1.4 for air
tc.gamma = gamma;

% Lax problem
lambda =10;

tc.rhoL = 0.445*lambda;
tc.uL   = 0.689*lambda ;
tc.pL   = 3.528*lambda;
tc.rhoR = 0.5*lambda;
tc.uR   = 0.0*lambda;
tc.pR   = 0.571*lambda;

tc.mL = tc.rhoL * tc.uL;
tc.EL = tc.pL / (tc.gamma - 1) + tc.rhoL * tc.uL^2 / 2;
tc.mR = tc.rhoR * tc.uR;
tc.ER = tc.pR / (tc.gamma - 1) + tc.rhoR * tc.uR^2 / 2;

% initial solution
tc.rho0 = @(x) tc.rhoL * (x < 0.0) + tc.rhoR * (x >= 0.0);
tc.u0   = @(x) tc.uL   * (x < 0.0) + tc.uR   * (x >= 0.0);
tc.p0   = @(x) tc.pL   * (x < 0.0) + tc.pR   * (x >= 0.0);
tc.m0   = @(x) tc.mL   * (x < 0.0) + tc.mR   * (x >= 0.0);
tc.E0   = @(x) tc.EL   * (x < 0.0) + tc.ER   * (x >= 0.0);

% constants required to compute p-star and exact solution
tc.aL = sqrt(tc.gamma * tc.pL ./ tc.rhoL);
tc.aR = sqrt(tc.gamma * tc.pR ./ tc.rhoR);
tc.AL = 2 / ((tc.gamma + 1) * tc.rhoL);
tc.BL = (tc.gamma - 1) / (tc.gamma + 1) * tc.pL;
tc.AR = 2 / ((tc.gamma + 1) * tc.rhoR);
tc.BR = (tc.gamma - 1) / (tc.gamma + 1) * tc.pR;

tc.TOL = 1e-6;

% Two-Rarefaction approximation
tc.pTR = ((tc.aL + tc.aR - (tc.gamma - 1) * (tc.uR - tc.uL) / 2) / ...
    (tc.aL / tc.pL^((tc.gamma - 1) / (2 * tc.gamma)) + ...
    tc.aR / tc.pR^((tc.gamma - 1) / (2 * tc.gamma))))^(2 * tc.gamma / (tc.gamma - 1));

% linearised solution
tc.pPV = (tc.pL + tc.pR) / 2 - (tc.uR - tc.uL) * (tc.rhoL + tc.rhoR) * (tc.aL + tc.aR) / 8;

% Two-Shock approximation
pHat = max(tc.TOL, tc.pPV);
gL = @(p)sqrt(tc.AL / (p + tc.BL));
gR = @(p)sqrt(tc.AR / (p + tc.BR));
tc.pTS = (gL(pHat) * tc.pL + gR(pHat) * tc.pR - (tc.uR - tc.uL)) / (gL(pHat) + gR(pHat));

% arithmetic mean
tc.pME = (tc.pL + tc.pR) / 2;

% domain of computation
tc.dm = [0, 1];

end
