function md = computeMeshData_line(msh)

isPeriodic = msh.bndTypes == 1;
if any(isPeriodic)
    IndFaceIDs = [msh.bndLFaces{isPeriodic}, msh.intLFaces];
else
    IndFaceIDs = msh.intLFaces;
end

md.nIntLFaces = length(IndFaceIDs);
md.intLFaces = cell(1, 3);
md.intLFaces{1, 1} = 2;
md.intLFaces{1, 2} = 1;
md.intLFaces{1, 3} = IndFaceIDs;

nBndTypes = length(msh.bndTypes);
md.nBndLFaces = msh.nBndLFaces;
md.bndLFaces = cell(2, nBndTypes);
for i = 1 : nBndTypes
    if (msh.bndTypes(i) == 1)
        md.nBndLFaces(i) = 0;
    else
        for lfn = 1 : 2
            md.bndLFaces{lfn, i} = msh.bndLFaces{i}(msh.faceNums(1, msh.bndLFaces{i}) == lfn);
        end
    end
end

end
