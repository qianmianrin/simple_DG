function [LU,label,dampingout] = applyprocess(U,dt,msh,md,bs,quad,tc,slo,pp,pro)

global limitertype


switch limitertype
    
    
    case{'jump_filter'}
         
        [LU,dampingout] = jumpfilter1DEuler(U,dt,msh,md,bs,quad,tc,pro);
        label = []; 
        LU = applyPPLimiter(LU, pp, tc, bs); 
    
      
    otherwise
        LU = U;
        label = [];
        dampingout = [];
end


end
