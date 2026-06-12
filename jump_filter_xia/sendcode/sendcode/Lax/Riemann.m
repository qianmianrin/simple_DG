%**************************************************************************
%                       problem description
%**************************************************************************
% Use EXRK-DG method for the one-dimensional Euler equations.
% pro = 6 is the only supported case: the Lax problem.
% rhoL = 0.445; uL = 0.698; pL = 3.528;
% rhoR = 0.5;   uR = 0;     pR = 0.571;
%
% basisType  : type of basis functions
% flux       : flux = 1 : Lax-Friedrichs flux
%            : flux = 2 : local Lax-Friedrichs flux
%            : flux = 3 : HLLC flux
% slo        : slo = 0 means no slope limiter
%              slo = 1 means TVD limiter
%              slo = 2 means TVB limiter
% pp         : pp = 0 means no positivity-preserving limiter
%              pp = 1 means positivity-preserving limiter
% pat        : pat = 0 means only plot at the final time
%            : pat = 1 means plot all the time
function [U0, t, pn, tc, quad1, bs,msh,dampingout] = ...
    Riemann(N, tp, pro, basisType, flux, slo, pp, pat)
%**************************************************************************
%                        some preparation work
%**************************************************************************
% tolerance for the matrix entries
mtol = 1.0e-10;

if (nargin < 3)
    error('Not enough arguments')
end

if pro ~= 6
    error('Only pro = 6 (Lax problem) is supported')
end

if (nargin < 4) || isempty(basisType)
    basisType = 102;
end
if basisType ~= 102
    error('Only basisType = 102 is supported')
end

if (nargin < 5) || isempty(flux)
    flux = 2;
end
if all(flux ~= 1 : 3)
    error('Wrong flux index')
end

if (nargin < 6) || isempty(slo)
    slo = 2;
end
if all(slo ~= 0 : 2)
    error('Wrong slope limiter index')
end

if (nargin < 7) || isempty(pp)
    pp = 0;
end
if all(pp ~= 0 : 1)
    error('Wrong positivity-preserving limiter index')
end

if (nargin < 8) || isempty(pat)
    pat = 0;
end
if all(pat ~= 0 : 1)
    error('Wrong pat index')
end

% Lax test case
tc = createTestCase(1.4, pro);

% quadrature rule for basisType = 102
k = 2;
quad1 = GaussQuadratureRule_line(3, 102);

% basis function set data
bs = setBasisFunctionSet_line(quad1, basisType);

% mesh
msh = setLineMesh_line(tc.dm, N, [2, 2], 101, 0, []);

% Augment mesh data
md = computeMeshData_line(msh);

% inverse of mass matrix at reference line
IME = inv(computeElementMatrix_refLine(0, 0, bs));
IME = mychop(IME, mtol);

% initial numerical solution
U0 = computeInitialSolution_line(msh, {tc.rho0, tc.m0, tc.E0}, bs, 1);
U0 = applyPPLimiter(U0, pp, tc, bs);

% CFL number
if (pp == 0)
    cfl = setCFLNumber(k);
else
    cfl = quad1.weights(1) / 2;
end

% points to plot
pn = setPlotPoints(msh, [], [], quad1);

%**************************************************************************
%                                 solve
%**************************************************************************
% integration in time using explicit Runge-Kutta methods
t = 0;
trouu = 11;
tstep = 0;
trouble = zeros(500, N);

Res = zeros(1, 1e+6);
Time = zeros(1, 1e+6);
nres = 0;
dampingout = [];

while(t < tp - 1.0e-12)
    dt = setdt(msh, U0, t, tc, quad1, bs, cfl, tp);
    [U, label, dampingout] = computeOneTimeStepEXRK(msh, md, U0, dt, flux, slo, pp, tc, bs, IME, quad1, pro);
    t = t + dt;
    tstep = tstep + 1;
    troub_label = find(label > 0);

    
    nres = nres + 1;
    Time(nres) = t;

    if mod(nres, 60000) == 0
        semilogy(1 : nres, Res(1 : nres), '-r', 'LineWidth', 2)
        str = strcat(int2str(nres), 'res');
        saveas(gcf, [str, '.fig']);
    end

    if trouu == 1
        if (~isempty(troub_label))
            trouble(tstep, 1:length(troub_label)) = troub_label;
        end
    end
 
    U0 = U;
end

 

 

end
