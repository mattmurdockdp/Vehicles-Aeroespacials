%% Problem 4: Compute the eigenfrequencies and eigenmodes of the model.
% Authors: Matheus and David
clear; 
load common_data.mat

n.modes = 100;

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

save constrained_modes.mat V_cN d_c inN dof2Restrict freq_u freq_c