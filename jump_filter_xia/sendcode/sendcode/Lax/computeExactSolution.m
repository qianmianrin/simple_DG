% exact solution for Riemann problem
function [rho, u, p] = computeExactSolution(x, t, tc)
    
if (t == 0)   
    rho = tc.rho0(x);
    u = tc.u0(x);
    p = tc.p0(x);   
    return   
end
p0 = max(tc.TOL, tc.pTS);
ps = computePStar(tc, p0);
us = computeUStar(tc, ps);

S = x / t;
if (ps > tc.pL)
    SL = tc.uL - tc.aL * sqrt((tc.gamma + 1) / (2 * tc.gamma) * ps / tc.pL + (tc.gamma - 1) / (2 * tc.gamma));
    rhosL = tc.rhoL * (ps / tc.pL + (tc.gamma - 1) / (tc.gamma + 1)) / ((tc.gamma - 1) / (tc.gamma + 1) * ps / tc.pL + 1);
else
    SHL = tc.uL - tc.aL;
    STL = us - tc.aL * (ps / tc.pL)^((tc.gamma - 1) / (2 * tc.gamma));
    rhosL = tc.rhoL * (ps / tc.pL)^(1 / tc.gamma);
    rhoLfan = tc.rhoL * (2 / (tc.gamma + 1) + (tc.gamma - 1) / (tc.gamma + 1) / tc.aL * (tc.uL - S)).^(2 / (tc.gamma - 1));
    uLfan = 2 / (tc.gamma + 1) * (tc.aL + (tc.gamma - 1) / 2 * tc.uL + S);
    pLfan = tc.pL * (2 / (tc.gamma + 1) + (tc.gamma - 1) / (tc.gamma + 1) / tc.aL * (tc.uL - S)).^(2 * tc.gamma / (tc.gamma - 1));
end

if (ps > tc.pR)
    SR = tc.uR + tc.aR * sqrt((tc.gamma + 1) / (2 * tc.gamma) * ps / tc.pR + (tc.gamma - 1) / (2 * tc.gamma));
    rhosR = tc.rhoR * (ps / tc.pR + (tc.gamma - 1) / (tc.gamma + 1)) / ((tc.gamma - 1) / (tc.gamma + 1) * ps / tc.pR + 1);
else
    SHR = tc.uR + tc.aR;
    STR = us + tc.aR * (ps / tc.pR)^((tc.gamma - 1) / (2 * tc.gamma));
    rhosR = tc.rhoR * (ps / tc.pR)^(1 / tc.gamma);
    rhoRfan = tc.rhoR * (2 / (tc.gamma + 1) - (tc.gamma - 1) / (tc.gamma + 1) / tc.aR * (tc.uR - S)).^(2 / (tc.gamma - 1));
    uRfan = 2 / (tc.gamma + 1) * (-tc.aR + (tc.gamma - 1) / 2 * tc.uR + S);
    pRfan = tc.pR * (2 / (tc.gamma + 1) - (tc.gamma - 1) / (tc.gamma + 1) / tc.aR * (tc.uR - S)).^(2 * tc.gamma / (tc.gamma - 1));  
end
    
if (ps > tc.pL)    
    if (ps > tc.pR)          
        rho = tc.rhoL * (S < SL) + rhosL * (S >= SL & S < us) +  rhosR * (S >= us & S < SR) + tc.rhoR * (S >= SR);
        u = tc.uL * (S < SL) + us * (S >= SL & S < SR) + tc.uR * (S >= SR);
        p = tc.pL * (S < SL) + ps * (S >= SL & S < SR) + tc.pR * (S >= SR);        
    else         
        rho = tc.rhoL * (S < SL) + rhosL * (S >= SL & S < us) + rhosR * (S >= us & S < STR) + rhoRfan .* (S >= STR & S < SHR) + tc.rhoR * (S >= SHR);
        u = tc.uL * (S < SL) + us * (S >= SL & S < STR) + uRfan .* (S >= STR & S < SHR) + tc.uR * (S >= SHR);
        p = tc.pL * (S < SL) + ps * (S >= SL & S < STR) + pRfan .* (S >= STR & S < SHR) + tc.pR * (S >= SHR);       
    end    
else
    if (ps > tc.pR)        
        rho = tc.rhoL * (S < SHL) + rhoLfan .* (S >= SHL & S < STL) + rhosL * (S >= STL & S < us) + rhosR * (S >= us & S < SR) + tc.rhoR * (S >= SR);
        u = tc.uL * (S < SHL) + uLfan .* (S >= SHL & S < STL) + us * (S >= STL & S < SR) + tc.uR * (S >= SR);
        p = tc.pL * (S < SHL) + pLfan .* (S >= SHL & S < STL) + ps * (S >= STL & S < SR) + tc.pR * (S >= SR);        
    else        
        rho = tc.rhoL * (S < SHL) + rhoLfan .* (S >= SHL & S < STL) + rhosL * (S >= STL & S < us) + rhosR * (S >= us & S < STR) + rhoRfan .* (S >= STR & S < SHR) + tc.rhoR * (S >= SHR);
        u = tc.uL * (S < SHL) + uLfan .* (S >= SHL & S < STL) + us * (S >= STL & S < STR) + uRfan .* (S >= STR & S < SHR) + tc.uR * (S >= SHR);
        p = tc.pL * (S < SHL) + pLfan .* (S >= SHL & S < STL) + ps * (S >= STL & S < STR) + pRfan .* (S >= STR & S < SHR) + tc.pR * (S >= SHR);        
    end   
end
                                 
end

%*********************************************************************************************
%                         subroutine to compute p-star                  
%*********************************************************************************************
function re = fL(tc, p)
if (p > tc.pL)
    re = (p - tc.pL) * sqrt(tc.AL / (p + tc.BL));
else
    re = 2 * tc.aL / (tc.gamma - 1) * ((p / tc.pL)^((tc.gamma - 1) / (2 * tc.gamma)) - 1);
end      
end

function re = fR(tc, p)
if (p > tc.pR)
    re = (p - tc.pR) * sqrt(tc.AR / (p + tc.BR));
else
    re = 2 * tc.aR / (tc.gamma - 1) * ((p / tc.pR)^((tc.gamma - 1) / (2 * tc.gamma)) - 1);
end      
end

function re = f(tc, p)
re = fL(tc, p) + fR(tc, p) + tc.uR - tc.uL;     
end

function re = dfL(tc, p)
if (p > tc.pL)
    re = sqrt(tc.AL / (p + tc.BL)) * (1 - (p - tc.pL) / (2 * (p + tc.BL)));
else
    re = (p / tc.pL)^(-(tc.gamma + 1) / (2 * tc.gamma)) / (tc.rhoL * tc.aL);
end      
end

function re = dfR(tc, p)
if (p > tc.pR)
    re = sqrt(tc.AR / (p + tc.BR)) * (1 - (p - tc.pR) / (2 * (p + tc.BR)));
else
    re = (p / tc.pR)^(-(tc.gamma + 1) / (2 * tc.gamma)) / (tc.rhoR * tc.aR);
end      
end

function re = df(tc, p)
re = dfL(tc, p) + dfR(tc, p);     
end

%copmpute p-star using Newton iteration
function [ps, iter] = computePStar(tc, p0)

p = p0 - f(tc, p0) / df(tc, p0);
p = max(tc.TOL, p); 
pm = (p0 + p) / 2;
CHA = abs(p - p0) / max(tc.TOL, pm);
iter = 1;

while (CHA >= tc.TOL)
    p0 = p;
    p = p0 - f(tc, p0) / df(tc, p0);
    p = max(tc.TOL, p); 
    pm = (p0 + p) / 2;
    CHA = abs(p - p0) / max(tc.TOL, pm);
    iter = iter + 1;    
end
ps = p;
end

%*********************************************************************************************
%                         subroutine to compute u-star                  
%*********************************************************************************************
function us = computeUStar(tc, ps)
us = (tc.uL + tc.uR) / 2 + (fR(tc, ps) - fL(tc, ps)) / 2;
end
