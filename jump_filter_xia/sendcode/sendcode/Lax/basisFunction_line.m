function r = basisFunction_line(x, type, index, der_order)

switch index
    case 1
        r = LegendrePolynomial(x, 0, der_order);
    case 2
        r = LegendrePolynomial(x, 1, der_order);
    case 3
        r = LegendrePolynomial(x, 2, der_order);
    otherwise
        error('Only three basis functions are available for basisType = 102')
end

end
