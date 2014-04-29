%% advect linearization
clear all;
clc;
close all;
addpath(sprintf('/Users/kevin/SkyDrive/KTH Work/LaTeX Reports',...
'/HW4-High Resolution shock-capturing methods/matlabfiles'));

%% Parameters
N = 200;
L = 10.;
dx = L/N; % Grid spacing
H = 1;
g = 9.8;
c = g; % Wave speed
tau = .5*dx/g; % Time Step
coeff = -tau/(2*dx);
nStep = 10*L/(c*tau);
i = 1:N;
%% Boundary Conditions
w = 0.1*L; % Width
x = (i-1/2)*dx; % Grid points
a = 1/5*H;
h = H + a*exp(-x.^2/(w^2));
m = zeros(1,N);
% Reflection Boundary Conditions
ip = i+1;
ip(N) = 1;
im = i-1;
im(1) = N;
%% Record the initial state
hplot(:,1) = h(:);
mplot(:,1) = m(:);
flux = (1/2*g*h.^2 + (m.^2)./h);
%% Loop
for iStep=1:(round(nStep))
    hOld=h;
    mOld=m;
    flux = (1/2*g*hOld.^2 + (mOld.^2)./hOld);
    h(1:N) = .5*(hOld(ip)+hOld(im)) + coeff*(mOld(ip)-mOld(im));
    m(1:N) = .5*(mOld(ip)+mOld(im)) + coeff*(flux(ip)-flux(im));
    hplot(:,iStep) = h(:);
    mplot(:,iStep) = m(:);
end

%% Plot
figure(1); clf;
plot(x,hplot(:,1),'-',x,h,'--');
legend('t=0  ','t=1');
xlabel('x');  ylabel('h(x,t)');
title(sprintf('$h(x,t)$ $\\Delta t=$%0.2g, $\\epsilon =$ %g',tau,a),...
    'Interpreter','latex')
%% h plot
figure(2)
for ip = 1:8:nStep
    clf;
    plot(x,hplot(:,1))
    hold on;
    plot(x,hplot(:,ip),'-');
    pause(.01)
end
%% m plot
% figure(3)
% for ip = 1:8:nStep
%     clf;
%     plot(x,mplot(:,1))
%     hold on;
%     plot(x,mplot(:,ip),'-');
%     pause(.01)
% end
%%
printYesNo = 0;
if printYesNo == 1
    saveFigurePath = '/Users/kevin/SkyDrive/KTH Work/LaTeX Reports/Stability HW3/Figures/';
    set(figure(1), 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('plot1%g',floor(tau^-1))]);
end