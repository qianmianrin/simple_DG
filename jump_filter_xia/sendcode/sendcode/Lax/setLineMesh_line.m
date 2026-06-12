function msh = setLineMesh_line(domain, N, bcs, type, maxLevel, funx)

if (nargin < 2)
    error('Not enough arguments')
end

if (length(domain) ~= 2) || (domain(2) <= domain(1))
    error('Wrong argument domain')
end

if (length(N) ~= 1)
    error('Wrong size of argument N')
end

if (nargin < 3) || isempty(bcs)
    bcs = ones(1, 2);
end
if (length(bcs) ~= 2)
    error('Wrong boundary conditions')
end

if (bcs(1) == 1)
    if (bcs(2) == 1)
        isPeriodicInX = true;
    else
        error('Wrong boundary condition in x direction')
    end
else
    if (bcs(2) == 1)
        error('Wrong boundary condition in x direction')
    else
        isPeriodicInX = false;
    end
end

if (nargin < 4) || isempty(type)
    type = 101;
end
if (type ~= 101) && (type ~= 102)
    error('wrong mesh type')
end

if (nargin < 5) || isempty(maxLevel) || (type == 101)
    maxLevel = 0;
end
if (type == 102) && (maxLevel == 0)
    maxLevel = 3;
end

if (nargin < 6) || isempty(funx)
    funx = @(x)x;
end
if abs(funx(0)) > 1.0e-12 || abs(funx(1) - 1) > 1.0e-12
    error('Wrong given scaling function in x direction')
end

xx = domain(1) + (domain(2) - domain(1)) * funx(linspace(0, 1, N + 1));
hx = diff(xx);

msh = getEmptyMesh;
msh.dm = domain;
msh.N = N;
msh.type = type;
msh.maxLevel = maxLevel;
msh.bndTypes = unique(bcs, 'stable');
msh.nElems = N;
msh.nFaces = N + 1;
if isPeriodicInX
    msh.nFaces = msh.nFaces - 1;
end
msh.nLElems = msh.nElems;
msh.nLFaces = msh.nFaces;
msh.LElems = 1 : msh.nElems;
msh.LFaces = 1 : msh.nFaces;

msh.elemCenter = xx(1 : N) + 0.5 * hx;
msh.elemLength = hx;
msh.elemSize = hx;
msh.elemFaces = zeros(2, msh.nElems);
msh.elemFaces(1, :) = 1 : N;
if isPeriodicInX
    msh.elemFaces(2, :) = [2 : N, 1];
else
    msh.elemFaces(2, :) = 2 : N + 1;
end
msh.elemJac = 0.5 * msh.elemSize;

msh.faceNormalx = ones(1, msh.nFaces);
msh.faceType = zeros(1, msh.nFaces);
msh.faceType(1) = bcs(1);
if ~isPeriodicInX
    msh.faceNormalx(1) = -1;
    msh.faceType(end) = bcs(2);
end

msh.faceElems = zeros(2, msh.nFaces);
msh.faceNums = zeros(2, msh.nFaces);
msh.faceElems(1, 2 : N) = 1 : N - 1;
msh.faceElems(2, 2 : N) = 2 : N;
msh.faceNums(1, 2 : N) = 2;
msh.faceNums(2, 2 : N) = 1;

if isPeriodicInX
    msh.faceElems(1, 1) = N;
    msh.faceElems(2, 1) = 1;
    msh.faceNums(1, 1) = 2;
    msh.faceNums(2, 1) = 1;
else
    msh.faceElems(1, 1) = 1;
    msh.faceNums(1, 1) = 1;
    msh.faceElems(1, N + 1) = N;
    msh.faceNums(1, N + 1) = 2;
end

msh.intLFaces = 2 : N;
msh.nIntLFaces = N - 1;
msh.bndLFaces = cell(1, length(msh.bndTypes));
msh.nBndLFaces = zeros(1, length(msh.bndTypes));
for i = 1 : length(msh.bndTypes)
    msh.bndLFaces{i} = find(msh.faceType == msh.bndTypes(i));
    msh.nBndLFaces(i) = length(msh.bndLFaces{i});
end

if (type == 102)
    msh.elemLevel = zeros(1, msh.nElems);
    msh.elemLID = 1 : msh.nElems;
    msh.elemParent = zeros(1, msh.nElems);
    msh.elemChildren = zeros(2, msh.nElems);
end

end
