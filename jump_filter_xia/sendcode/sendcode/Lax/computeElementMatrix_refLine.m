function EM = computeElementMatrix_refLine(trial_der, test_der, bs, mtol)

if (nargin < 1) || isempty(trial_der)
    trial_der = 0;
end

if (nargin < 2) || isempty(test_der)
    test_der = 0;
end

if (trial_der ~= 0) && (trial_der ~= 1)
    error('Not implemented derivative order of trial basis function')
end

if (test_der ~= 0) && (test_der ~= 1)
    error('Not implemented derivative order of test basis function')
end

if (nargin < 3) || isempty(bs) || ~isstruct(bs)
    error('basisType = 102 requires a basis function set struct')
end

if (nargin < 4) || isempty(mtol)
    mtol = 1.0e-12;
end

if ~strcmpi(bs.refGeom, 'line')
    error('Wrong reference geometry for basis functions to evaluate on')
end
if bs.type ~= 102
    error('Only basisType = 102 is supported')
end

if (bs.elemPointsType == 102) && (2 * bs.nep - 3 >= 2 * bs.deg)
    EM = bs.phitw{test_der + 1} * bs.phi{trial_der + 1};
else
    quad = GaussQuadratureRule_line(bs.deg + 2, 102);
    phi_trial = basisFunctionSet_line(quad.points, bs.type, trial_der);
    phi_test = basisFunctionSet_line(quad.points, bs.type, test_der);
    phitw_test = (quad.weights .* phi_test)';
    EM = phitw_test * phi_trial;
end

EM = mychop(EM, mtol);

end
