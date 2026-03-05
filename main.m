%% Problema 0: Preprocessing

close all;
clear;
load model_data_clean.mat

% Struct de mides 'n'
n.nodes = length(node_coords);  % Nombre total de nodes
n.dims  = 6;                    % Nombre de DOF per cada node 
n.dof   = n.dims*n.nodes;       % Nombde total de DOF

%% Problema 1: càlcul estructural estàtic considerant només els 3 suports fixes

% Indexadors dels nodes de les diferents condicions de contorn
[inD, ...   %   Dirichlet
 inN] ...   %   Neumann
= findBCIndicies(case_control_sets.subcase_0_SET_2, n);

% Forces
g_node = [0;0;-9.81;0;0;0];
gVec   = repmat(g_node,n.nodes,1);
FT     = MAAX * gVec;

% Massa total
m = 0;
for i = 1:n.nodes
    idx = (i-1)*6 + 1;
    m   = m + MAAX(idx,idx);
end

FN    = FT(inN);  % Vector de forces amb condicions de Neumann   
FDext = FT(inD);  % Vector de forces externes per trobar les reaccions de Dirichlet

% Sistema d'equacions
uD  = zeros(length(inD),1);
KNN = KAAX(inN,inN);
KND = KAAX(inN,inD);
KDD = KAAX(inD,inD);
KDN = KAAX(inD,inN);
uN  = KNN\(FN-KND*uD);  
FD  = KDD*uD+KDN*uN;

% Reaccions 
rF = FD-FDext;

% Comparar reaccions amb el pes
FGravity = m*9.81;
fprintf('\n---- Comparació de reaccions amb el pes ----\n');
fprintf('Pes total (N):            %.6f\n', FGravity);
fprintf('Pes/3 (per suport) (N):   %.6f\n\n', FGravity/3);

fprintf('Reacció suport 1 (N):     %.6f\n', rF(3));
fprintf('Reacció suport 2 (N):     %.6f\n', rF(9));
fprintf('Reacció suport 3 (N):     %.6f\n', rF(15));

fprintf('\nSuma reaccions (N):       %.6f\n', rF(3)+rF(9)+rF(15));