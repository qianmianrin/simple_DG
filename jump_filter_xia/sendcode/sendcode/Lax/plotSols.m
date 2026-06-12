function plotSols(U, t, pn, tc, quad, bs,msh,dampingout)

% exact values to plot
[exact_rho,~,~] = computeExactSolution(pn.exa, t, tc);
 
% numerical values to plot
numerical_rho = 0.5 * reshape(quad.weights' * (bs.phi{1} * U(:, 3 * pn.numElems - 2)), [length(pn.num), 1]);
 
% subplot(2, 2, 1)
figure(1)
plot(pn.exa, (exact_rho), '-b', pn.num, (numerical_rho), '-r');
xlabel('x', 'FontSize', 16);
ylabel('density', 'FontSize', 16);
legend('\rho', '\rho_h', 'FontSize', 16);
title('average density');

end

