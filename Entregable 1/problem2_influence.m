%% Problem 2 – Influence matrix
% Authors: Matheus and David 
clear; 
load common_data.mat

% BCs 
dof2Restrict = 1:3;
[inD, inN] = findBCIndicies(case_control_sets.subcase_0_SET_2, n, dof2Restrict);

% Polinomis Zernike
n.noll = 100;
Z = zeros(n.nodes,n.noll);
for j=1:n.noll
    [~,~,~,Z(:,j),~] = zernike_noll(node_coords(:,1)/radi, node_coords(:,2)/radi, j);
end

% Calcul dels coeficients d'influència
n.act = length(case_control_sets.subcase_0_SET_1);
U = zeros(n.nodes,n.act);
C = zeros(n.noll,n.act);

for i = 1:n.act
    FT = zeros(n.dof,1);
    dof = node2DOF(case_control_sets.subcase_0_SET_1(i),n,1:3);
    FT(dof(3)) = 1;

    [Uaux,~,~] = solveStatics(inD,inN,FT,KAAX,n);
    U(:,i) = Uaux(3:6:end);

    C(:,i) = (Z \ U(:,i)) * 1000; % um/N
end

% Resultats per l'actuador 13 pels modes de 4 a 10
fprintf("\n---- Coeficients de la matriu d'influencia per l'actuador 13, dels modes 4 al 10----\n");
for k = 4:10
    fprintf('Mode %d : %d µm/N\n',k,C(k,13));
end

save influence_matrix.mat C Z U n