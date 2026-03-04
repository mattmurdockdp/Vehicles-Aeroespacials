function DOFList = node2DOF(nodeList, n)
    % Descripció:
    %   Funció que converteix una llista de nodes en una llista que conté
    %   tots els seus graus de llibertat, en notació de Voight
    %
    % Arguments:
    %   nodeList    [nDOF_subset, 1]
    %   n           [struct]            Struct de mides, interessa n.dims
    %
    % Sortides:
    %   DOFList     [nDIMS*nDOF_subset, 1] = 
    %       [n_1^1, n_1^2, ... n_{nDOF_subset}^nDIMS]
    %
    % Última modificació:
    %   04.03.2026
    %
    % Funció verificada!

    % nodeList = (1:3)';
    DOFList = zeros(n.dims*length(nodeList), 1);

    for i = 1:length(DOFList)
        subIdx = floor((i-1)/n.dims)+1;      % Índex pel subset de nodes
        % disp(subIdx);                      % Comprovar que subIdx funciona
        j = mod(i-1,n.dims)+1;               % Índex de graus de llibertat
        % disp(j);                           % Comprovar que j funciona
        
        DOFList(i) = j + (nodeList(subIdx)-1)*n.dims;
    end    

end