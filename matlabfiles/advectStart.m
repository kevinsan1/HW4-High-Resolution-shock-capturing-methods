%% advect linearization
clear all;
clc;
close all;
addpath(['/Users/kevin/SkyDrive/KTH Work/LaTeX Reports',...
'/HW4-High Resolution shock-capturing methods/matlabfiles/']);

%% Parameters
N = 10;
L = 1.;
dx = L/N; % Grid spacing
H = 1;
g = 9.8;
c = g; % Wave speed
tau = .5*dx/g; % Time Step
coeff = -tau/(2*dx);
nStep = L/(c*tau);
i = 1:N;
ghostCellOneSide = 2;
%% Boundary Conditions
w = 0.1*L; % Width
x = (i-1/2)*dx; % Grid points
a = 1/5*H;
h = H + a*exp(-(x-L/2).^2/(w^2));
m = zeros(1,N + ghostCellOneSide*2);
% Boundary Conditions
h = padarray(h',ghostCellOneSide)';
%% BC 
% needs to be changed if ghostCellOneSide is changed
h(1) = h(4);
h(2) = h(3);
m(1) = m(4)*(-1);
m(2) = m(3)*(-1);
h(N+2) = h(N-1);
h(N+1) = h(N);
m(N+2) = m(N-1)*(-1);
m(N+1) = m(N)*(-1);
ip = 
%% Record the initial state
hplot(:,1) = h(:);
mplot(:,1) = m(:);
flux = (1/2*g*h.^2 + (m.^2)./h);
%% Loo
for iStep=1:(round(nStep))
    hOld=h;
    mOld=m;
    flux = (1/2*g*hOld.^2 + (mOld.^2)./hOld);
    h(1:N) = .5*(hOld(ip)+hOld(i)) + coeff*(mOld(ip)-mOld(i));
    m(1:N) = .5*(mOld(ip)+mOld(i)) + coeff*(flux(ip)-flux(i));
    hplot(:,iStep) = h(:);
    mplot(:,iStep) = m(:);
end

%% Plot
figure(1); clf;
plot(x,hplot(:,1),'-',x,h,'--');
legend('t=0  ','t=1');
xlabel('x');  ylabel('h(x,t)');
title(sprintf('$h(x,t)$ $\\Delta t=$%0.2g,',...
    ' $\\epsilon =$ %g',tau,a),...
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
    saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/'...
        'LaTeX Reports/HW4-High Resolution shock-capturing'...
        'methods/Figures/'];
    set(figure(1), 'PaperPositionMode', 'auto');
    print('-depsc2', [saveFigurePath ...
        sprintf('plot1%g',floor(tau^-1))]);
end