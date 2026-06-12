function bs = setBasisFunctionSet_line(quad, type)

% reference geometry to evaluate
bs.refGeom = 'line';

% type of basis functions
bs.type = type;

% the maximum polynomial degree w.r.t. one variable
bs.deg = mod(type, 10);

% points type for evaluation
bs.elemPointsType = quad.type;

% Evalute basis functions and their derivatives at Gauss nodes of the
% reference domain [-1, 1]
bs.phi = {basisFunctionSet_line(quad.points, type, 0), ...
          basisFunctionSet_line(quad.points, type, 1)};

% Evalute basis functions at face of the refence domain [-1, 1]
bs.phi_face = {basisFunctionSet_line(-1, type, 0), ...
               basisFunctionSet_line(1, type, 0)}; 
bs.phi_facederx = {basisFunctionSet_line(-1, type, 1), ...
               basisFunctionSet_line(1, type, 1)}; 
bs.phi_facederxx = {basisFunctionSet_line(-1, type, 2), ...
               basisFunctionSet_line(1, type, 2)}; 
bs.phi_facederxxx = {basisFunctionSet_line(-1, type, 3), ...
               basisFunctionSet_line(1, type, 3)};            
% Weights multiply basis functions and their derivatives componentwisely to 
% prepare for integration
bs.phitw = {(quad.weights .* bs.phi{1})', ...
            (quad.weights .* bs.phi{2})'};

% transpose of basis function evaluations at face
bs.phitw_face = {bs.phi_face{1}', bs.phi_face{2}'};

% number of basis functions and its square
bs.nb = size(bs.phi{1}, 2);
bs.nb2 = bs.nb * bs.nb;

% number of Gauss points at element
bs.nep = quad.np;

end

