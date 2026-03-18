%% MAIN – Problem A
% Authors: Matheus and David
% Aquest script executa cada mòdul de forma individual. S'utilitzen archius
% del tipus .mat amb les variables necessàries guardades en cada mòdul.
% Es pot comentar cada mòdul en funció del que es vulgui veure.

close all; clear; clc;

% Problem 0 – Preprocessing
problem0_preprocess;

% Problem 1 – Gravity sag
fprintf('\n')
disp("------ Problem 1: Gravity Sag ------")
problem1_gravity;

% Problem 2 – Influence matrix
fprintf('\n')
disp("------ Problem 2: Influence Matrix ------")
problem2_influence;

% Problem 3 – Force optimization
fprintf('\n')
disp("------ Problem 3: Force Optimization ------")
problem3_optimization;

% Problem 4 – Modal analysis
fprintf('\n')
disp("------ Problem 4: Modal Analysis ------")
problem4_modes;

% Problem 5 – Dynamic compliance
fprintf('\n')
disp("------ Problem 5: Dynamic Compliance ------")
disp("Veure Plot")
problem5_dynamic;
  