function r = LegendrePolynomial(x, poly_deg, der_order)

if (poly_deg < 0) || (poly_deg > 2)
    error('Only polynomial degrees 0, 1, and 2 are supported')
end

switch der_order
    case 0
        switch poly_deg
            case 0
                r = ones(size(x));
            case 1
                r = x;
            case 2
                r = x.^2 - 1 / 3;
        end
    case 1
        switch poly_deg
            case 0
                r = zeros(size(x));
            case 1
                r = ones(size(x));
            case 2
                r = 2 * x;
        end
    case 2
        switch poly_deg
            case {0, 1}
                r = zeros(size(x));
            case 2
                r = 2 * ones(size(x));
        end
    case 3
        r = zeros(size(x));
    otherwise
        error('Only derivative orders 0, 1, 2, and 3 are supported')
end

end
