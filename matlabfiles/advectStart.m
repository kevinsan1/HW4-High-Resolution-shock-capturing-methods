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
c = g*sqrt(H); % Wave speed
tau = .5*dx/g; % Time Step
coeff = -c*tau/(2*dx);
nStep = 2*L/(c*tau);
nCells = 1:N;
ghostCellOneSide = 2;
cellsPlusGhost = 1:(N+2*ghostCellOneSide);
%% Boundary Conditions
w = 0.1*L; % Width
x = (nCells-1/2)*dx; % Grid points
a = 1/5*H;
h = H + a*exp(-(x-L/2).^2/(w^2));
h = padarray(h',ghostCellOneSide)';
% Boundary Conditions
m = g*sqrt(h)*0;
%% BC
% needs to be changed if ghostCellOneSide is changed
h(1) = h(4);
h(2) = h(3);
m(1) = m(4)*(-1);
m(2) = m(3)*(-1);
h(N+4) = h(N+1);
h(N+3) = h(N+2);
m(N+4) = m(N+1)*(-1);
m(N+3) = m(N+2)*(-1);
%%
iWithGhost = ((1+ghostCellOneSide):(N+ghostCellOneSide));
ip = iWithGhost+1;
im = iWithGhost-1;
%% Record the initial state
hplot(:,1) = h(:);
mplot(:,1) = m(:);
flux = (1/2*g*h.^2 + (m.^2)./h);
%% Loo
for iStep=1:(round(nStep))
    hOld=h;
    mOld=m;
    flux = (1/2*g*hOld.^2 + (mOld.^2)./hOld);
    h(3:(N+2)) = .5*(hOld(ip)+hOld(im)) + coeff*(mOld(ip)-mOld(iWithGhost));
    m(3:(N+2)) = .5*(mOld(ip)+mOld(im)) + coeff*(flux(ip)-flux(iWithGhost));
    h(1) = h(4);
    h(2) = h(3);
    m(1) = m(4)*(-1);
    m(2) = m(3)*(-1);
    h(N+4) = h(N+1);
    h(N+3) = h(N+2);
    m(N+4) = m(N+1)*(-1);
    m(N+3) = m(N+2)*(-1);
    hplot(:,iStep) = h(:);
    mplot(:,iStep) = m(:);
end

%% Plot
% figure(1); clf;
% plot(x,hplot(3:(N+2),1),'-',x,h,'--');
% legend('t=0  ','t=1');
% xlabel('x');  ylabel('h(x,t)');
% title(sprintf('$h(x,t)$ $\\Delta t=$%0.2g,',...
%     ' $\\epsilon =$ %g',tau,a),...
%     'Interpreter','latex')
%% h plot1
figure(2)
for ip = 1:2:nStep
    clf;
    plot(x,hplot((3:(N+2)),1))
    hold on;
    plot(x,hplot((3:(N+2)),ip),'-');
    pause(.1)
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