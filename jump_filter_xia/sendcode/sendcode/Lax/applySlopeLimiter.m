% lim = 0 :  no limiter
% lim = 1 : TVD limiter
% lim = 2 : TVB limiter
function U = applySlopeLimiter(msh, md, U, slo, tc, bs)

switch slo
    
    % no limiter
    case 0
        
        return 
                   
    % slope limiter(TVD or TVB)
    case {1, 2} 
        
        if bs.type ~= 102
            error('Only basisType = 102 is supported')
        end
        
        % Set up a small number
        eps = 1e-7; 
                        
        % Compute the difference of mean values in neighboring elements
        um   = reshape(U(1, :), [3, msh.nLElems]);
        fdum = zeros(3, msh.nLElems);
        bdum = zeros(3, msh.nLElems);
        
        faceIDs = md.intLFaces{1, 3};
        leLIDs  = msh.faceElems(1, faceIDs);
        reLIDs  = msh.faceElems(2, faceIDs);
        fdum(:, leLIDs) = um(:, reLIDs) - um(:, leLIDs);
        bdum(:, reLIDs) = fdum(:, leLIDs); 
        
%         faceID = md.bndLFaces{1, 1};  
%         leLID = msh.faceElems(1, faceID);
%         bdum(:, leLID) = um(:, leLID) - [tc.rhoL; tc.mL; tc.EL].*100 ; 
% 
%         faceID = md.bndLFaces{2, 1};
%         leLID = msh.faceElems(1, faceID);
%         fdum(:, leLID) = [tc.rhoR; tc.mR; tc.ER].*100  - um(:, leLID); 
                 
        % Compute the modified slop of P1 part, 
        % i.e. the modified second coefficient
        ux     = reshape(U(2, :), [3, msh.nLElems]);
        ux_mod = zeros(3, msh.nLElems);
        if (slo == 1)
            for ie = 1 : msh.nLElems        
                % Compute the eigenmatrix and its inverse 
                pm = computePressure(um(1, ie), um(2, ie), um(3, ie), tc);
                cm = computeSpeedOfSound(um(1, ie), pm, tc);
                Hm = computeEnthalpy(um(1, ie), um(3, ie), pm);
                Rm = computeEigenmatrix(um(2, ie) / um(1, ie), cm, Hm);
                RIm = computeEigenmatrixInv(um(2, ie) / um(1, ie), cm, tc);
                
                ux_mod(:, ie) = Rm * minmod(RIm * ux(:, ie), RIm * fdum(:, ie), RIm * bdum(:, ie));
            end  
        else
            M = 100;
            h = msh.elemLength(:, msh.LElems);
            for ie = 1 : msh.nLElems        
                % Compute the eigenmatrix and its inverse 
                pm = computePressure(um(1, ie), um(2, ie), um(3, ie), tc);
                cm = computeSpeedOfSound(um(1, ie), pm, tc);
                Hm = computeEnthalpy(um(1, ie), um(3, ie), pm);
                Rm = computeEigenmatrix(um(2, ie) / um(1, ie), cm, Hm);
                RIm = computeEigenmatrixInv(um(2, ie) / um(1, ie), cm, tc);
                
                ux_mod(:, ie) = Rm * minmod_bar(RIm * ux(:, ie), RIm * fdum(:, ie), RIm * bdum(:, ie), M, h(1, ie));
            end                       
        end
        ux_mod = ux_mod(:)';
        
        % Reconstruction
        isModified = abs(ux_mod - U(2, :)) >= eps;
        U(2, isModified) = ux_mod(isModified);
        if (bs.deg ~= 1)           
            U(3 : end, isModified) = 0;            
        end
        
        return
end
                                                           
end
