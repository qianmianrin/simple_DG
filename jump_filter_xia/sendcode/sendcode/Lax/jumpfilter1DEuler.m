function [LU,dampingout] = jumpfilter1DEuler(U,dt,msh,md,bs,quad,tc,pro)

% written by Lei Wei 2023.11.18

h = msh.elemLength;

rho = U(:,1:3:end);
rhou = U(:,2:3:end);
Ener = U(:,3:3:end);
polydegree = size(rho,1) - 1;

% radius
rhohave = repmat(0.5 * (quad.weights' * (bs.phi{1} * rho)), size(rho,1), 1);
rhouhave = repmat(0.5 * (quad.weights' * (bs.phi{1} * rhou)), size(rho,1), 1);
uhave = rhouhave ./ rhohave;
Ehave = repmat(0.5 * (quad.weights' * (bs.phi{1} * Ener)), size(rho,1), 1);
pave = computePressure(rhohave, rhouhave, Ehave, tc);
cave = computeSpeedOfSound(rhohave, pave, tc);
betaj = abs(uhave) + abs(cave);

dampingrho = calculate_jump(rho, md, msh, bs, quad, 'rho', pro);
dampingrhou = calculate_jump(rhou, md, msh, bs, quad, 'rhou', pro);
dampingEner = calculate_jump(Ener, md, msh, bs, quad, 'Ener', pro);

enthapy = (Ehave(1,:) + pave(1,:)) ./ rhohave(1,:);

para = 1.0  ;
cf = para * betaj(1,:);
cf = cf ./ enthapy ./ 2 ;

damping = zeros(size(rho));
damping(2,:) = max([dampingrho(2,:); dampingrhou(2,:); dampingEner(2,:)]);

if polydegree == 2
    damping(3,:) = max([dampingrho(3,:); dampingrhou(3,:); dampingEner(3,:)]);
end

if polydegree == 3
    damping(3,:) = max([dampingrho(3,:); dampingrhou(3,:); dampingEner(3,:)]);
    damping(4,:) = max([dampingrho(4,:); dampingrhou(4,:); dampingEner(4,:)]);
end

h = repmat(h, size(rho,1), 1);
damping = damping .* cf;
damping = -damping * dt ./ h;

Lrho = exp(damping) .* rho;
Lrhou = exp(damping) .* rhou;
LEner = exp(damping) .* Ener;

LU = U;
LU(:,1:3:end) = Lrho;
LU(:,2:3:end) = Lrhou;
LU(:,3:3:end) = LEner;

dampingout = -damping ./ dt;

end
