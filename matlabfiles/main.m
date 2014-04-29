%% HW 4
%
clear all;
clc;
close all;
addpath(sprintf('/Users/kevin/SkyDrive/KTH Work/LaTeX Reports',...
'/HW4-High Resolution shock-capturing methods/matlabfiles'));
%% Parameters
L = 1; % Total length
N = 80; % Number of grid spaces
dx = L/N; % Grid spacing
H = 1;
a = H/5;
w = 0.1*L;
g = 9.8;
c = g; % Wave speed
tau = 0.1;
grid_numbers = 1:N; 
% grid_numbers = [1,2,3,4,5] if N = 5
x = (grid_numbers - 1/2)*dx; % values of x
% x(1) = (1 - 1/2) * dx = (   0.5   ) * L/N = L/(2N)
% x(N) = (N - 1/2) * dx = ( N - 1/2 ) * L/N = L - L/(2N)
%% Initialize h and m
h = H + a*exp(-(x-L/2).^2/(w^2));
m = zeros(1,N);
%% Main loop
for i = 1:
%% Plot h and m
% Plot h
figure(1)
plot(h)
% Plot m
figure(2)
plot(m)