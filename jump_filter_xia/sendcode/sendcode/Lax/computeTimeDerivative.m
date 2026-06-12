% Compute residual
function Ut = computeTimeDerivative(msh, md, U, flux, tc, bs, IME, pro)

if pro ~= 6
    error('Only pro = 6 (Lax problem) is supported')
end

% Evaluate physical quantities and flux functions
ul = reshape(bs.phi_face{1} * U, [3, msh.nLElems]);
ur = reshape(bs.phi_face{2} * U, [3, msh.nLElems]);
pl = computePressure(ul(1, :), ul(2, :), ul(3, :), tc);
pr = computePressure(ur(1, :), ur(2, :), ur(3, :), tc);
Fl = computeF(ul(1, :), ul(2, :), ul(3, :), pl);
Fr = computeF(ur(1, :), ur(2, :), ur(3, :), pr);

%**************************************************************************
%                      element contributions
%**************************************************************************
if (bs.type == 100)
    Ut = zeros(1, 3 * msh.nLElems);
else
    u = bs.phi{1} * U;
    p = computePressure(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), tc);
    F = computeF(u(:, 1 : 3 : end), u(:, 2 : 3 : end), u(:, 3 : 3 : end), p);
    Ut = bs.phitw{2} * reshape(F, [bs.nep, 3 * msh.nLElems]);
end

%**************************************************************************
%                       internal face contributions
%**************************************************************************
% dissipation coefficient for the LF flux
alpha = 1;
if (flux == 1)
    alpha = computeDissCoe;
end

faceIDs = md.intLFaces{1, 3};
leLIDs  = msh.faceElems(1, faceIDs);
reLIDs  = msh.faceElems(2, faceIDs);
nf      = length(faceIDs);

F_hat = computeInviscidFlux(ur(:, leLIDs), ul(:, reLIDs), pr(leLIDs), pl(reLIDs), Fr(:, leLIDs), Fl(:, reLIDs), ones(1, nf), flux, tc);
Ut(:, (-2 : 0)' + 3 * leLIDs) = Ut(:, (-2 : 0)' + 3 * leLIDs) - bs.phitw_face{2} * reshape(F_hat, [1, 3 * nf]);
Ut(:, (-2 : 0)' + 3 * reLIDs) = Ut(:, (-2 : 0)' + 3 * reLIDs) + bs.phitw_face{1} * reshape(F_hat, [1, 3 * nf]);

%**************************************************************************
%                       boundary face contributions
%**************************************************************************
% left
faceID = md.bndLFaces{1, 1};
leLID  = msh.faceElems(1, faceID);
u_ext = [tc.rhoL; tc.mL; tc.EL];
F_hat = computeInviscidFlux(ul(:, leLID), u_ext, pl(leLID), pl(leLID), -Fl(:, leLID), -Fl(:, leLID), -1, flux, tc);
Ut(:, (-2 : 0)' + 3 * leLID) = Ut(:, (-2 : 0)' + 3 * leLID) - bs.phitw_face{1, 1} * F_hat';

% right
faceID = md.bndLFaces{2, 1};
leLID  = msh.faceElems(1, faceID);
u_ext = [tc.rhoR; tc.mR; tc.ER];
F_hat = computeInviscidFlux(ur(:, leLID), u_ext, pr(leLID), pr(leLID), Fr(:, leLID), Fr(:, leLID), 1, flux, tc);
Ut(:, (-2 : 0)' + 3 * leLID) = Ut(:, (-2 : 0)' + 3 * leLID) - bs.phitw_face{1, 2} * F_hat';

%**************************************************************************
%                        post-processing
%**************************************************************************
% Take care of the mass matrix
Ut = (IME * Ut) ./ repelem(msh.elemJac(:, msh.LElems), 1, 3);

%**************************************************************************
%                          subroutine
%**************************************************************************
% Compute the dissipation coefficient for the LF flux
    function alpha = computeDissCoe
        alpha = max(max(computeEigenMax(ul(1, :), ul(2, :), pl, tc), computeEigenMax(ur(1, :), ur(2, :), pr, tc)));
    end

% Compute inviscid numerical flux
    function F_hat = computeInviscidFlux(UL, UR, pl, pr, FL, FR, n, flux, tc)
        
        switch flux
            case 1
                F_hat = computeLFFlux(UL, UR, FL, FR, alpha);
            case 2
                F_hat = computeLLFFlux(UL, UR, pl, pr, FL, FR, tc);
            case 3
                F_hat = computeHLLCFlux(UL, UR, pl, pr, FL, FR, n, tc);
        end
        
    end

end
