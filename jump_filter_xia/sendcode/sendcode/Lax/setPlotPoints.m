function pn = setPlotPoints(msh, stride, np, quad)

if isempty(stride)
    if isempty(np)
        v = 1 : msh.nLElems;
    else
        if (np > msh.nLElems)
            np = msh.nLElems;
        end        
        v = round(linspace(1, msh.nLElems, np));
    end
else
    v = 1 : stride : msh.nLElems;
    if v(end) ~= msh.nLElems
        v = [v, msh.nLElems];
    end     
end

% Sort the each leaf element according to its center
[sct, I] = sort(msh.elemCenter(1, msh.LElems));
LElems = msh.LElems(I);

% points to plot at for exact solutions
pn.exa = reshape(sct + 0.5 * msh.elemLength(:, msh.LElems) .* quad.points, [1, msh.nLElems * quad.np]);

% points to plot at for numerical solutions
pn.num      = sct(v);
pn.numElems = LElems(v);

end