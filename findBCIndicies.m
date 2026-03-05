function [inD, inN] = findBCIndicies(subset, n, nod2Restrict)
    % Descripció:
    %   Funció que genera els vectors amb els índexs de cada tipus de condició
    %   de contorn.
    %
    % Arguments:
    %   subset  [nNodes_subset, 1]  Llista de nodes amb condició de
    %           Dirichlet
    %   n       [struct]            Struct de mides
    %
    % Sortides:
    %   inD     [nDOF_D, 1] Índexs dels graus de llibertat amb condició de
    %           Dirichlet
    %   inN     [nDOF_N, 1] Índexs dels graus de llibertat amb condició de
    %           Neumann
    %
    % Última modificació:
    %   05.03.2026
    %
    % Funció verificada!

    inD = node2DOF(subset,n,nod2Restrict);
    inT = 1:n.dof;
    inN = setdiff(inT, inD);

end