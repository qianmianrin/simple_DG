function quad = GaussQuadratureRule_line(np, type)



quad.type = type;
quad.np = np;

switch np
    case 3
        quad.points = [-1; 0; 1];
        quad.weights = [1 / 3; 4 / 3; 1 / 3];
    case 4
        quad.points = [-1; -sqrt(1 / 5); sqrt(1 / 5); 1];
        quad.weights = [1 / 6; 5 / 6; 5 / 6; 1 / 6];
    otherwise
        error('Only 3- and 4-point Gauss-Lobatto rules are needed for basisType = 102')
end

end
