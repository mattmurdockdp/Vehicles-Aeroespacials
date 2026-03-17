function b = reshapeForMETA(a,nDims)
% Descripció:
    %   Funció que passa un vector de [n.dofs,1] -> [n.nodes,n.dims] perquè
    %   el format encaixi amb la funció donada fillhdf.m
    %
    % Arguments:
    %   a       [n.dofs,1]          Vector genèric
    %   nDim    [escalar]           # graus de llibertat per node
    %
    % Sortides:
    %   b       [n.nodes,n.dims]    Vector genèric amb format canviat
    %
    % Última modificació:
    %   16.03.2026
    % Funció verificada!
    
    % Es calculen les dimensions específiques pel vector entrat, així també
    % serveix per vectors amb nodes restringits (menor nombre de GDL)
    nDof = length(a);
    nNodes = nDof/nDims;
    
    b = zeros(nNodes,nDims);
    for j=1:nDims
        b(:,j) = a(j:nDims:nDof);
    end

end