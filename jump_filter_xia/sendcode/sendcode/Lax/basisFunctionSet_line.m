function r = basisFunctionSet_line(x, type, der_order)


x = x(:);
r = [basisFunction_line(x, type, 1, der_order), ...
     basisFunction_line(x, type, 2, der_order), ...
     basisFunction_line(x, type, 3, der_order)];

end
