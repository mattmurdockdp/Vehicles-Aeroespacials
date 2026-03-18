%% Problem 5: Dynamic compliance of actuator 13 and usable band-width
% Authors: Matheus and David
clear;
load common_data.mat
load constrained_modes.mat   % V_cN, d_c, inN, dof2Restrict

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
f     = linspace(0,2000,500);    % Hz
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
xlim([0,1000])
legend('A(f)','+3 dB','Location','best');
set(gcf,'Color','w'); saveas(gcf,'dynamic_compliance.png');