%% Problema 3 - Force optimization
% Authors: Matheus and David 
% a) Calcular el force vector amb el mètode dels minims quadrats
clear; 
load common_data.mat
load influence_matrix.mat   % C, Z, U, n

% Definim Zt tenint en compte els valors de la taula (en nm RMS) i que no optimitzarem
% els noll index < 4
Zt = zeros(n.noll,1);

Zt(4)  = 1300; Zt(5)  = -650; Zt(6)  = 650;
Zt(7)  = -350; Zt(8)  = 350;  Zt(9)  = -180;
Zt(10) = 180;  Zt(11) = -120; Zt(12) = 70;
Zt(13) = -70;  Zt(14) = -45;  Zt(15) = 45;

% Definim només amb noll index > 4
C_opt  = C(4:100, :).*1000; % nm/N rms
Zt_opt = Zt(4:100);         % nm rms

% Resolució per mètode de mínims quadrats
f = C_opt \ Zt_opt; % [n.act x 1]

% Calculem els coeficients de zernike després de compensació
Zc = C_opt*f - Zt_opt;

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

save optimization_results.mat f Zc Zt_opt C_opt sigma_total sigma_100 sigma_hf
