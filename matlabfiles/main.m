%% HW 4
%
clear all;
clc
addpath('/Users/kevin/SkyDrive/KTH Work/LaTeX Reports/HW4-High Resolution shock-capturing methods/matlabfiles');
%% Parameters
L = 10; % Total length
N = 10; % Number of grid spaces
dx = L/N; % Grid spacing
grid_numbers = 1:N; 
% grid_numbers = [1,2,3,4,5] if N = 5
x = (grid_numbers - 1/2)*dx; % values of x
% x(1) = (1 - 1/2) * dx = (   0.5   ) * L/N = L/(2N)
% x(N) = (N - 1/2) * dx = ( N - 1/2 ) * L/N = L - L/(2N)
