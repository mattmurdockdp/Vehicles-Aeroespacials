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
dof2Restrict = 1:3;

[inD, ...   %   Dirichlet
    inN] ...   %   Neumann
    = findBCIndicies(case_control_sets.subcase_0_SET_2, n, dof2Restrict);

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
FT = zeros(n.dof, 1);
C = zeros(n.noll,n.act);

for i = 1:n.act
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
    C(:,i) = (Z\W)*1000;    % mm/N
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
C_opt  = C(4:100, :)*1000; % nm/N rms
Zt_opt = Zt(4:100);       % nm rms

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
n.modes = 200;

% a) Unconstrained.
[V_u, d_u] = eigs(KAAX, MAAX, n.modes, 'smallestabs');
freq_u     = sqrt(diag(d_u))/(2*pi);
% La freqüència surt en números complexos.

disp('-----------------------------');
disp('M O D E S   P R O P I S');
disp(' No restringit');

for i = 1:8
    ratio = abs(freq_u(i+1)/freq_u(i));
    fprintf('Ratio de freq. modes %d i %d: %d Hz \n', i+1, i, ratio);
end
% NOTA: revisar si el ràtio ha de ser gran o petit, de cara al report


% b) Constrained.
% Restringim tots els graus de llibertat dels tres suports
dof2Restrict = 1:6;
[inD, ...   %   Dirichlet
    inN] ...   %   Neumann
    = findBCIndicies(case_control_sets.subcase_0_SET_2, n, dof2Restrict);

% Seleccionem només els graus de llibertat de Neumann, és a dir, els no
% restringits.
[V_cN, d_c] = eigs(KAAX(inN,inN), MAAX(inN,inN), n.modes, 'smallestabs');
freq_c      = sqrt(diag(d_c))/(2*pi);

% Els graus de llibertat de Dirichlet tindran desplaçament nul. Cal
% considerar-los pel format de META.
V_c        = zeros(n.dof,n.modes);
V_c(inN,:) = V_cN;

% Exporting the 6th mode for visualization in META
V_c6th = reshapeForMETA(V_c(:,6),6);
fillhdf('h5template.h5', '6thModeConstrained.h5', V_c6th*1e-3); % Passem el vector propi en [m]

disp(' Suports com a nodes restringits');

for i = 1:8
    ratio = abs(freq_c(i+1)/freq_c(i));
    fprintf('Ratio de freq. modes %d i %d: %d Hz \n', i+1, i, ratio);
end

%% Problem 5: Dynamic compliance of actuator 13 and usable band-width

zeta = 0.002; % modal damping ratio

% Constrained modes
Vmodal = V_cN;          % Aquí no cal considerar els nodes de Dirichlet
Dmodal = d_c;
Mmodal = MAAX(inN,inN);

% freqüències modals
omega_i = real(sqrt(abs(diag(Dmodal))));   % rad/s
f_i = omega_i/(2*pi);                      % Hz
nModes = size(Vmodal,2);

% Es calculen els paràmetres modals
m_i = zeros(nModes,1);
k_i = zeros(nModes,1);
b_i = zeros(nModes,1);

for i=1:nModes
    phi    = Vmodal(:,i);
    m_i(i) = phi' * Mmodal * phi.*1000; % passem a kg
    k_i(i) = m_i(i) .* (omega_i(i).^2);
    b_i(i) = 2 .* m_i(i) .* omega_i(i) .* zeta;
end

% Es busca el dof per la força vertical de l'actuador 13
node_act13  = case_control_sets.subcase_0_SET_1(13);
dof_act13   = node2DOF(node_act13,n,dof2Restrict);
dof_act13_z = dof_act13(3);

% Terme amplitud (tenint en compte que la força és 1 N)
phi_i = Vmodal(dof_act13_z,:).';

% Fem el mapa de freqüències, el qual ha de ser de 0 a 1000 hz, però per
% estudiar-se s'ha de fer fins al doble, és a dir, fins a 2000 hz.
f     = linspace(0,2000,3000);    % Hz
omega = 2*pi*f;                % rad/s
H     = zeros(size(omega));        % dynamic compliance

for j = 1:length(omega)
    w     = omega(j);
    denom = - (w^2) .* m_i + k_i + 1i .* (b_i .* w);   % (nModes x 1)
    xi_r  = phi_i ./ denom;
    H(j)  = sum(phi_i .* xi_r)*1000; % formula 36 de la teoria, convertim a micras
end

% static compliance
H0 = sum( (phi_i.^2) ./ k_i )*1000;  % amb w = 0

% Guany
A = 20*log10( abs(H) ./ abs(H0) );

% Es calcula el usable bandwidth com el bandwidth pel cual A <= 3 dB
idx_ok   = find(A <= 3);
last_idx = max(idx_ok); % volem la freq maxima
f_bw     = f(last_idx);

% --- print results ---
fprintf('\nDynamic compliance at actuator %d (z DOF):\n', 13);
fprintf('Static compliance H(0) = %.6e um/N\n', H0);
fprintf('Usable bandwidth (A <= +3 dB): f_bw = %.3f Hz\n\n', f_bw);

% --- plots ---
figure('Name','Dynamic compliance and amplification','NumberTitle','off','Units','normalized','Position',[0.1 0.1 0.6 0.6]);
subplot(2,1,1);
semilogx(f, abs(H),'b','LineWidth',1.2);
xlabel('Frequency (Hz)'); ylabel('|H| (um/N)');
title(sprintf('Dynamic compliance at actuator %d (z DOF)', 13));
grid on;

subplot(2,1,2);
plot(f, A,'k','LineWidth',1.0);
hold on;
yline(3,'r--','+3 dB','LineWidth',1);
xlabel('Frequency (Hz)'); ylabel('Amplification A(f) (dB)');
xlim([0 2000]);
ylim([min(A)-1 max(A)+1]);
grid on;
legend('A(f)','+3 dB','Location','best');



