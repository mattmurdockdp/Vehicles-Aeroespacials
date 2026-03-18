%% Problem 1 – Gravity sag
% Authors: Matheus and David 
load common_data.mat

dof2Restrict = 1:3;
[inD, inN] = findBCIndicies(case_control_sets.subcase_0_SET_2, n, dof2Restrict);

g_node = [0;0;-g;0;0;0];
FT     = repmat(g_node,n.nodes,1);
FT     = MAAX * FT;

% Massa total
m = 0;
for i = 1:n.nodes
    idx = (i-1)*6 + 1;
    m   = m + MAAX(idx,idx);
end

[U, F, rF] = solveStatics(inD, inN, FT, KAAX, n);
U = reshapeForMETA(U,6);
fillhdf('h5template.h5','gravitySagU.h5',U*1e-3);

% Comparar reaccions amb el pes
FGravity = m*g;
fprintf('\n---- Comparació de reaccions amb el pes ----\n');
fprintf('Pes total (N):            %.6f\n', FGravity);
fprintf('Pes/3 (per suport) (N):   %.6f\n\n', FGravity/3);

fprintf('Reacció suport 1 (N):     %.6f\n', rF(3));
fprintf('Reacció suport 2 (N):     %.6f\n', rF(6));
fprintf('Reacció suport 3 (N):     %.6f\n', rF(9));

fprintf('\nSuma reaccions (N):       %.6f\n', rF(3)+rF(6)+rF(9));


save gravity_results.mat U rF inD inN