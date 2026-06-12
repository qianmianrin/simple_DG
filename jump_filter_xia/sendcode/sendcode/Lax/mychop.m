function M = mychop(M, eps)

if (nargin < 2) || isempty(eps)
    eps = 1.0e-10;
end

M = M .* (M < -eps | M > eps);

end
