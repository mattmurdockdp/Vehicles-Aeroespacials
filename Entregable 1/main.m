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

% Indexadors dels nodes de les diferents condicions de contorn
dof2Restrict = 1:3;

% DOFs per les BC
[inD, ...   %   Dirichlet
    inN] ...   %   Neumann
    = findBCIndicies(case_control_sets.subcase_0_SET_2, n, dof2Restrict);

%% Problema 1: càlcul estructural estàtic considerant només els 3 suports fixes

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

% Solve system
[U, F, rF] = solveStatics(inD, inN, FT, KAAX, n);
U = reshapeForMETA(U,6);
fillhdf('h5template.h5', 'gravitySagU.h5', U*1e-3); % Passem desplaçaments en [m]

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
C = zeros(n.noll,n.act);

for i = 1:n.act
    % Inicialitzem per cada actuador per no acumular forçes d'actuadors que
    % no toquen
    FT = zeros(n.dof, 1);
    % Trobem els dofs d'aquest actuador
    actuatorDof = node2DOF(case_control_sets.subcase_0_SET_1(i), n, dof2Restrict);
    % Nota: nod2Restrict aquí es pot utilitzar el que es vulgui, sempre
    % i quan sigui robust amb el que s'indexa dins de FT aquí a sota
    FT(actuatorDof(3)) = 1; % 1 N positiu en l'eix z

    % idxZ = repmat([0; 0; 1; 0; 0; 0], n.nodes, 1);
    [Uaux, F, rF] = solveStatics(inD, inN, FT, KAAX, n);

    % Guardem components Z
    U(:,i) = Uaux(3:6:end); % mm/N

    % Apartat c, calcular la matriu d'influència
    W      = U(:,i);
    C(:,i) = (Z\W)*1000;    % um/N
end

% Resultats per l'actuador 13 pels modes de 4 a 10
fprintf("\n---- Coeficients de la matriu d'influencia per l'actuador 13, dels modes 4 al 10----\n");
for k = 4:10
    fprintf('Mode %d : %d µm/N\n',k,C(k,13));
end


%% Problema 3

% a) Calcular el force vector amb el mètode dels minims quadrats

% Definim Zt tenint en compte els valors de la taula (en nm RMS) i que no optimitzarem
% els noll index < 4

Zt = zeros(n.noll,1);

Zt(4)  = 1300; Zt(5)  = -650; Zt(6)  = 650;
Zt(7)  = -350; Zt(8)  = 350;  Zt(9)  = -180;
Zt(10) = 180;  Zt(11) = -120; Zt(12) = 70;
Zt(13) = -70;  Zt(14) = -45;  Zt(15) = 45;

% Definim només amb noll index > 4
C_opt  = C(4:100, :).*1000; % nm/N rms
Zt_opt = Zt(4:100);        % nm rms

% Resolució per mètode de mínims quadrats
f = C_opt \ Zt_opt; % [n.act x 1]

% Calculem els coeficients de zernike després de compensació
Zc = C_opt*f-Zt_opt;

% i) force vector for actuator from 13 to 17
fprintf("\n---- Vector de forces pels actuadors del 13 al 17 ----\n");
for k = 13:17
    fprintf('Actuador %d : %d N\n',k,f(k));
end
fprintf("\n")

% ii) The total compensation residual RMS, rejecting the rigid body modes

% calculem el desplaçament nodal
Wc = U*f;  % n.nodes x 1 en mm

% s'ha de tenir en compte que l'efecte dels rigid body (noll<=3) s'ha d'eliminar
Z_rb     = Z(:,1:3);        % Z en mm
coeff_rb = Z_rb \ Wc;       % C pels modes rb en mm/N
Wc_rigid = Z_rb * coeff_rb; % Desplaçaments dels modes rb en mm
Wc_noRB  = Wc - Wc_rigid;   % Eliminem desplaçaments dels modes rb
Wc_noRB  = Wc_noRB*1e6;     % Passem a nm

% RMS total
sigma_total = sqrt(sum(Wc_noRB.^2));
fprintf('Error total rms (sense modes rigid body)  : %.6f nm\n',sigma_total);

% iii) Compensation residual captured by the 100 modes not considering rb
% modes
sigma_100 = sqrt(sum(Zc.^2));
fprintf('Residual capturat pels modes 4..100 de Zernike: %.6f nm\n', sigma_100);

% iv) High frequency (n>100) residual
sigma_hf = sqrt(sum(Wc_noRB.^2)-sum(Zc.^2));
fprintf('Residual de freqüència alta (RMS Noll>100) = %.6f nm\n', sigma_hf);

%% Problem 4: Compute the eigenfrequencies and eigenmodes of the model.

n.modes = 30;

% ---------------------------------------------------------------- %
% a) Unconstrained.
[V_u, d_u] = eigs(KAAX, MAAX, n.modes, 'smallestabs');
freq_u     = sqrt(diag(d_u))/(2*pi);   % Hz

disp('-----------------------------');
disp('M O D E S   P R O P I S');
disp(' No restringit');

% Print first 8 frequencies (Hz)
fprintf('\nPrimeres freqüències (unconstrained):\n');
for i = 1:8
    fprintf('f_%d = %.3f Hz\n', i, freq_u(i));
end

% Print ratios (como en el main)
fprintf('\nRatios de freqüència (unconstrained):\n');
for i = 1:8
    ratio = abs(freq_u(i+1)/freq_u(i));
    fprintf('Ratio modes %d i %d: %.3f\n', i+1, i, ratio);
end

% ---------------------------------------------------------------- %
% b) Constrained.
% Restringim tots els graus de llibertat dels tres suports
dof2Restrict = 1:6;
[inD, inN] = findBCIndicies(case_control_sets.subcase_0_SET_2, n, dof2Restrict);

% Seleccionem només els graus de llibertat de Neumann
[V_cN, d_c] = eigs(KAAX(inN,inN), MAAX(inN,inN), n.modes, 'smallestabs');
freq_c      = sqrt(diag(d_c))/(2*pi);   % Hz

% Els graus de llibertat de Dirichlet tindran desplaçament nul
V_c        = zeros(n.dof,n.modes);
V_c(inN,:) = V_cN;

% Exporting the 6th mode for visualization in META
V_c6th = reshapeForMETA(V_c(:,6),6);
fillhdf('h5template.h5', '6thModeConstrained.h5', V_c6th*1e-3);

disp(' ');
disp(' Suports com a nodes restringits');

% Print first 8 constrained frequencies
fprintf('\nPrimeres freqüències (constrained):\n');
for i = 1:8
    fprintf('f_%d = %.3f Hz\n', i, freq_c(i));
end

% Print ratios (como en el main)
fprintf('\nRatios de freqüència (constrained):\n');
for i = 1:8
    ratio = abs(freq_c(i+1)/freq_c(i));
    fprintf('Ratio modes %d i %d: %.3f\n', i+1, i, ratio);
end

%% Problem 5: Dynamic compliance of actuator 13 and usable band-width
zeta = 0.002; % modal damping ratio

% Constrained modes
Vmodal = V_cN;
Dmodal = d_c;
Mnn    = MAAX(inN,inN).*1000; % (t -> kg)
Knn    = KAAX(inN,inN).*1000; % (N/mm -> N/m)

% freqüències modals
omega_i = real(sqrt(abs(diag(Dmodal))));   % rad/s
nModes  = size(Vmodal,2);

% Fem el mapa de freqüències, el qual ha de ser de 0 a 1000 hz, però per
% estudiar-se s'ha de fer fins al doble, és a dir, fins a 2000 hz.
f     = linspace(0,2000,100);    % Hz
w     = 2*pi*f;                  % rad/s

% Inicialitzem els paràmetres modals i les forces unitàries per l'actuador 13
phi                        = Vmodal;
force_vector               = zeros(length(inN), 1);
node_act13                 = case_control_sets.subcase_0_SET_1(13);  % actuador 13
dof_act13                  = node2DOF(node_act13, n, dof2Restrict);  
dofz_global                = dof_act13(3);                           % DOF en z
idx                        = find(inN == dofz_global);               % posicio del DOF
force_vector(idx)          = 1;                                      % força unitaria aplicada al DOF (1 N)
dynamics_displacements     = zeros(length(f), nModes);
static_displacements       = zeros(length(f), nModes);

for i = 1:nModes
    % Massa modal (kg)
    m_i = phi(:,i)' * Mnn * phi(:,i);
    % Rigidesa modal (N/m)
    k_i = phi(:,i)' * Knn * phi(:,i);
    % Damper modal
    b_i = 2 * zeta * omega_i(i) * m_i;
    % Força modal
    f_modal = phi(:,i)' * force_vector;
    % Dynamic compliance (m)
    dynamics_displacements(:,i) = f_modal ./(k_i - m_i*w.^2 + 1i*b_i.*w);
    % Static compliance (m)
    static_displacements(:,i) = f_modal / k_i; 
end

% Nomes projectem a l'index de l'actuador ja que es un sistema desacoplat,
% aixi MATLAB no peta
H  = phi(idx, :) * dynamics_displacements.';   % 1 × nFreq
H0 = phi(idx, :) * static_displacements.';     % 1 × nFreq

% Guany
A = 20*log10( abs(H) ./ abs(H0) );

% --- plot ---
figure(1);
plot(f, A,'k','LineWidth',1.0);
hold on;
yline(3,'r--','+3 dB','LineWidth',1);
xlabel('Frequency (Hz)'); ylabel('Amplification A(f) (dB)');
grid on;
legend('A(f)','+3 dB','Location','best');



