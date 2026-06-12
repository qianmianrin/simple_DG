function [damping] = calculate_jump(u,md,msh,bs,~,~,~)
% u represents some conserved variable, e.g. rho, rhou, Ener.
polydegree = size(u,1) - 1;

faceIDs = md.intLFaces{1, 3};
leLIDs = msh.faceElems(1, faceIDs);
reLIDs = msh.faceElems(2, faceIDs);
h = msh.elemLength;

leLIDs = [1 leLIDs];
reLIDs = [1 reLIDs];

ul = bs.phi_face{1} * u;
ur = bs.phi_face{2} * u;

leftvertexjump  = abs(ur(leLIDs) - ul(reLIDs));
rightvertexjump = [leftvertexjump(2:end) leftvertexjump(1)];
[leftvertexjump,rightvertexjump] = modifyboundary(leftvertexjump,rightvertexjump);
jumpmoment1 = leftvertexjump + rightvertexjump;

ulder = 2 ./ h .* (bs.phi_facederx{1} * u);
urder = 2 ./ h .* (bs.phi_facederx{2} * u);

leftvertexjumpder  = abs(urder(leLIDs) - ulder(reLIDs));
rightvertexjumpder = [leftvertexjumpder(2:end) leftvertexjumpder(1)];
[leftvertexjumpder,rightvertexjumpder] = modifyboundary(leftvertexjumpder,rightvertexjumpder);
jumpmoment = jumpmoment1 + (leftvertexjumpder + rightvertexjumpder) .* 2 .* h;

damping = zeros(size(u));
damping(2,:) = jumpmoment;

if polydegree == 2
    ulderder = (2 ./ h).^2 .* (bs.phi_facederxx{1} * u);
    urderder = (2 ./ h).^2 .* (bs.phi_facederxx{2} * u);
    
    leftvertexjumpderder  = abs(urderder(leLIDs) - ulderder(reLIDs));
    rightvertexjumpderder = [leftvertexjumpderder(2:end) leftvertexjumpderder(1)];
    [leftvertexjumpderder,rightvertexjumpderder] = modifyboundary(leftvertexjumpderder,rightvertexjumpderder);

    jumpmoment = jumpmoment + (leftvertexjumpderder + rightvertexjumpderder) .* 6 .* h.^2;
    damping(3,:) = jumpmoment;
end

if polydegree == 3
    ulderder = (2 ./ h).^2 .* (bs.phi_facederxx{1} * u);
    urderder = (2 ./ h).^2 .* (bs.phi_facederxx{2} * u);
    
    leftvertexjumpderder  = abs(urderder(leLIDs) - ulderder(reLIDs));
    rightvertexjumpderder = [leftvertexjumpderder(2:end) leftvertexjumpderder(1)];
    [leftvertexjumpderder,rightvertexjumpderder] = modifyboundary(leftvertexjumpderder,rightvertexjumpderder);

    jumpmoment = jumpmoment + (leftvertexjumpderder + rightvertexjumpderder) .* 6 .* h.^2;
    damping(3,:) = jumpmoment;
    
    ulderderder = (2 ./ h).^3 .* (bs.phi_facederxxx{1} * u);
    urderderder = (2 ./ h).^3 .* (bs.phi_facederxxx{2} * u);
    
    leftvertexjumpderderder  = abs(urderderder(leLIDs) - ulderderder(reLIDs));
    rightvertexjumpderderder = [leftvertexjumpderderder(2:end) leftvertexjumpderderder(1)];
    [leftvertexjumpderderder,rightvertexjumpderderder] = modifyboundary(leftvertexjumpderderder,rightvertexjumpderderder);

    jumpmoment = jumpmoment + (leftvertexjumpderderder + rightvertexjumpderderder) .* 12 .* h.^3;
    damping(4,:) = jumpmoment;
end

end

function [a,b] = modifyboundary(a,b)

a(:,1) = 0;
b(:,end) = 0;

end
