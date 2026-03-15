%% Problema 0: Preprocessing

close all;
clear;
load model_data_clean.mat

% Constants
g = 9.81*1e3;               % [mm/s^2]
radi = 75;                  % [mm]

% Struct de mides 'n'
n.nodes = length(node_coords);  % Nombre total de nodes
n.dims  = 6;                    % Nombre de DOF per cada node 
n.dof   = n.dims*n.nodes;       % Nombde total de DOF

%% Problema 1: càlcul estructural estàtic considerant només els 3 suports fixes

% Indexadors dels nodes de les diferents condicions de contorn
nod2Restrict = 1:3;

[inD, ...   %   Dirichlet
 inN] ...   %   Neumann
= findBCIndicies(case_control_sets.subcase_0_SET_2, n, nod2Restrict);

% Forces
g_node = [0;0;-g;0;0;0];
gVec   = repmat(g_node,n.nodes,1);
FT     = MAAX * gVec;

% Massa total
m = 0;
for i = 1:n.nodes
    idx = (i-1)*6 + 1;
    m   = m + MAAX(idx,idx);
end

[U, F, rF] = solveStatics(inD, inN, FT, KAAX, n);

% Comparar reaccions amb el pes
FGravity = m*g;
fprintf('\n---- Comparació de reaccions amb el pes ----\n');
fprintf('Pes total (N):            %.6f\n', FGravity);
fprintf('Pes/3 (per suport) (N):   %.6f\n\n', FGravity/3);

fprintf('Reacció suport 1 (N):     %.6f\n', rF(3));
fprintf('Reacció suport 2 (N):     %.6f\n', rF(6));
fprintf('Reacció suport 3 (N):     %.6f\n', rF(9));

fprintf('\nSuma reaccions (N):       %.6f\n', rF(3)+rF(6)+rF(9));

%% Problema 2

% Apartat a, desplaçament dels nodes segons el Polinomi de Zernike
n.noll = 100;
Z = zeros(n.nodes, n.noll);
for j=1:n.noll              
    [~, ~, ~, Z(:, j), ~] = zernike_noll(node_coords(:,1)/radi, node_coords(:,2)/radi, j);
end

% Apartat b, efecte aïllat de cadascun dels actuadors
n.act = length(case_control_sets.subcase_0_SET_1);
U = zeros(n.nodes, n.act);
FT = zeros(n.dof, 1);

for i = 1:n.act
    actuatorDof = node2DOF(case_control_sets.subcase_0_SET_1(i), n, nod2Restrict);
        % Nota: nod2Restrict aquí es pot utilitzar el que es vulgui, sempre
        % i quan sigui robust amb el que s'indexa dins de FT aquí a sota
    FT(actuatorDof(3)) = 1; % 1 N positiu en l'eix z

    % idxZ = repmat([0; 0; 1; 0; 0; 0], n.nodes, 1);
    [Uaux, F, rF] = solveStatics(inD, inN, FT, KAAX, n);

    % Guardem components Z
    U(:,i) = Uaux(3:6:end);

end

% Apartat c, calcular la matriu d'influència

C = zeros(n.noll,n.act);

for i=1:n.act
    W = U(:,i);
    C(:,i) = (Z\W)*1000;
end

% Resultats per l'actuador 13 pels modes de 4 a 10
fprintf("\n---- Coeficients de la matriu d'influencia per l'actuador 13, dels modes 4 al 10----\n");
for k = 4:10
    fprintf('Mode %d : %d µm/N\n',k,C(k,13));
end


%% Problema 3

% a) Calcular el force vector amb el mètode dels minims quadrats

% Definim Zt tenint en compte els valors de la taula i que no optimitzarem
% els noll index < 4

Zt = zeros(n.noll,1);

Zt(4) = 1300; Zt(5)  = -650; Zt(6)  = 650;
Zt(7) = -350; Zt(8)  = 350;  Zt(9)  = -180;
Zt(10) = 180; Zt(11) = -120; Zt(12) = 70;
Zt(13) = -70; Zt(14) = -45; Zt(15)  = 45;

% Definim només amb noll index > 4  

C_opt = C(4:100, :)*1000; % nm/N rms
Zt_opt = Zt(4:100);       % nm rms

% Resolució per mètodes quadrats
f = C_opt \ Zt_opt;

% Calculem els coeficients de zernike després de compensació
Zc = C_opt*f-Zt_opt;

% i) force vector for actuator from 13 to 17
fprintf("\n---- Vector de forces pels actuadors del 13 al 17----\n");
for k = 13:17
    fprintf('Actuador %d : %d N\n',k,f(k));
end
fprintf("\n")

% ii) The total compensation residual RMS, rejecting the rigid body modes
sigma_total = sqrt(sum((Zc(:)).^2));
fprintf('Error total rms : %d N\n \n',sigma_total);

% iii) Compensation residual captured by the 100 modes

% iv) High frequency (n>100) residual

%% Problem 4 Compute the eigenfrequencies and eigenmodes of the model.

% a) Unconstrained.


% b) Constrained.



