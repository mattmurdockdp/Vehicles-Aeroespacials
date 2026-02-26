%% Problema 1

load model_data_clean.mat
support_node_list = case_control_sets.subcase_0_SET_2;

% Definim la matriu de fixnodes i calculem l'array d'index de Dirichlet
fixnodes = zeros(3*length(support_node_list),3);
inD = zeros(length(fixnodes),1);
i = 1;

for inode = 1:length(support_node_list)
    for dof=1:3
        index = dof + (inode-1)*3; 
        fixnodes(index,:) = [
            support_node_list(inode),dof,0
            ]; 
        % Array of Dirichlet dofs
        inD(i) = dof + (fixnodes(index,1)-1)*3;
        i = i+1;
    end
end

% Dirichlet displacements
uD = fixnodes(:,3);

% Array d'index de Neumann
iT = 1:length(node_coords)*3;
iT = iT';








