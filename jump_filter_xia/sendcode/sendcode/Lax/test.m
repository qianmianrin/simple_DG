clear; clc
format long
Globals1D

%******************parameters setting**************************************
% number of elements
N = 400;
% type of basis functions
basisType = 102;
% limiter type
limitertype = 'jump_filter';
% flux type, 1, 2 or 3
flux = 2;
% Lax problem
pro = 6;
% time to print
tp = 0.038;

% slope limiter type
slo = 1 + 1;

% use positivity preserving limiter or not
pp = 1;
% plot all the time or not
pat = 0;
%**************************************************************************
[U0, t, pn, tc, quad1, bs, msh, dampingout] = Riemann(N, tp, pro, basisType, flux, slo, pp, pat);

plotSols(U0, t, pn, tc, quad1, bs, msh);
