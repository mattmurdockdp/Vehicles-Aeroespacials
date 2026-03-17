%% Problem 0 – Preprocessing
% Authors: Matheus and David 
clear; clc;

load model_data_clean.mat

g    = 9.81*1e3;   % mm/s^2
radi = 75;         % mm

n.nodes = length(node_coords);
n.dims  = 6;
n.dof   = n.nodes*n.dims;

save common_data.mat g radi n node_coords KAAX MAAX case_control_sets