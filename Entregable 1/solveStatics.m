function [U, F, rF] = solveStatics(inD, inN, FT, KAAX, n)
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

    % Ensamblar
    U = zeros(n.dof, 1);
        U(inD) = uD;
        U(inN) = uN;
    F = zeros(n.dof, 1);
        F(inD) = FD;
        F(inN) = FN;
end