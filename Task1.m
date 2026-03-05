%% Problema 1
clear;
load model_data_clean.mat
support_node_list = case_control_sets.subcase_0_SET_2;

% Definim la matriu de fixnodes i calculem l'array d'index de Dirichlet
fixnodes = zeros(3*length(support_node_list),3);
iD = zeros(length(fixnodes),1);
i = 1;

for inode = 1:length(support_node_list)
    for dof=1:3
        index = dof + (inode-1)*3; 
        fixnodes(index,:) = [
            support_node_list(inode),dof,0
            ]; 
        % Array of Dirichlet dofs
        iD(i) = dof + (fixnodes(index,1)-1)*3;
        i = i+1;
    end
end

% Dirichlet displacements
uD = fixnodes(:,3);

% Array d'index de Neumann
iT = 1:length(node_coords)*3;
iT = iT';

iN = setdiff(iT, iD);


% Vector de forces global
Nnodes = length(node_coords);
Ndof = 3*Nnodes;
gVec = repmat([0;0;-9.81], Nnodes, 1); % Vector de gravetat
mVec = zeros(Ndof,1);

for i = 1:Nnodes
    for j = 1:3
        idx = (i-1)*6 + 1;
        mVec(3*(i-1)+j) = MAAX(idx,idx);
    end
end

FT = gVec.*mVec;

% Vector de forces dels nodes amb condició de Neumann
FN = FT(iN);

FDext = FT(iD);

% Sistema d'equacions
uD  = zeros(length(iD),1);
KNN = KAAX(iN,iN);
KND = KAAX(iN,iD);
KDD = KAAX(iD,iD);
KDN = KAAX(iD,iN);
uN  = KNN\(FN-KND*uD);  
FD  = KDD*uD+KDN*uN;

% Reaccions 
rF = FD-FDext;
