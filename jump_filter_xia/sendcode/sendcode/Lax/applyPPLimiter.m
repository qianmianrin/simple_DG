% pp = 0 : no limiter
% pp = 1 : pp limiter
function U = applyPPLimiter(U, pp, tc, bs)

switch pp
    
    % no limiter
    case 0
        
        return                  
       
    % PP limiter (positivity-preserving)
    case 1
        
        if bs.type ~= 102
            error('Only basisType = 102 is supported')
        end
        
        % Evaluate mean values of pressure
        pm = computePressure(U(1, 1 : 3 : end), U(1, 2 : 3 : end), U(1, 3 : 3 : end), tc); 

        % Set up a small number
        eps = min(1.0e-9, min(min(U(1, 1 : 3 : end), pm)));

        % Evaluate density before modification at Gauss-Lobatto points
        rho = bs.phi{1} * U(:, 1 : 3 : end);

        % Compute the scaling parameters theta1 for desity
        rho_min = min(rho);
        theta1 = min((U(1, 1 : 3 : end) - eps) ./ (U(1, 1 : 3 : end) - rho_min), 1);

        % U after scaling desity
        U(2 : end, 1 : 3 : end) = theta1 .* U(2 : end, 1 : 3 : end);

        % Evaluate values after modifying density
        u = bs.phi{1} * U;
        p = computePressure(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), tc);        

        % Compute the scaling parameters theta2 for pressure 
        % For points where pressure are less than eps, we obtain t_alpha
        % by solving quadratic equations
        t_alpha = ones(size(p));
        pos     = find(p < eps);
        col     = ceil(pos / bs.nep);
        if (~isempty(pos))
            rhom_pos = U(1, 3 * col - 2)';
            mm_pos   = U(1, 3 * col - 1)';
            Em_pos   = U(1, 3 * col)';
            pm_pos   = pm(col)';
            rho_pos  = getMatEntries(u(:, 1 : 3 : end), pos);
            m_pos    = getMatEntries(u(:, 2 : 3 : end), pos);
            E_pos    = getMatEntries(u(:, 3 : 3 : end), pos);
            p_pos    = p(pos);
            
            AA = (tc.gamma - 1) * (mm_pos .* m_pos - rhom_pos .* E_pos - rho_pos .* Em_pos);
            aa = rhom_pos .* pm_pos + rho_pos .* p_pos + AA;
            bb = -AA - 2 * rhom_pos .* pm_pos - eps * (rho_pos - rhom_pos);
            cc = rhom_pos .* (pm_pos - eps);
            t_alpha(pos) = (-bb - sqrt(bb.^2 - 4 * aa .* cc)) ./ (2 * aa);
            
            theta2 = min(t_alpha);
            % final U 
            U(2 : end, :) = repelem(theta2, 1, 3) .* U(2 : end, :); 
        end

        return 
end

function entries = getMatEntries(mat, pos)
    entries = mat(pos);
end

end
