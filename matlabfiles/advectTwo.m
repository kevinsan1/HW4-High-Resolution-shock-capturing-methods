%% advect linearization
clear all;

%% Parameters
N = 101;
L = 10.;
dx = L/N;               % Grid spacing
H = 1;
g = 9.8;
c = sqrt(g*H);        % Wave speed
tau = .2*dx/c;             % Time Step
coeff = -c*tau/(2*dx);
nStep = L/(c*tau);
i = 1:N;
%% Boundary Conditions
w = 0.4;              % Width
x = (i-1/2)*dx;    % Grid points
h = H + .1*exp(-(x-L/2).^2/(w^2));
m = h*c*tau;
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
for iStep=1:nStep
    hOld=h;
    mOld=m;
    flux = (1/2*g*hOld.^2 + (mOld.^2)./hOld);
    h(1:N) = .5*(hOld(ip)+hOld(im)) + 1/c*coeff*(mOld(ip)-mOld(im));
    m(1:N) = .5*(mOld(ip)+mOld(im)) + c*coeff*(flux(ip)-flux(im));
    hplot(:,iStep) = h(:);
    mplot(:,iStep) = m(:);
end

%% Plot
figure(1); clf;
plot(x,hplot(:,1),'-',x,h,'--');
legend('t=0  ','t=1');
xlabel('x');  ylabel('h(x,t)');
title(sprintf('$h(x,t)$'),...
    'Interpreter','latex')
%%
figure(2)

for ip = 1:4:nStep
    plot(x,hplot(:,1))
    hold on;
    plot(x,hplot(:,ip),'-');
    pause(.3)
    clf
end
%%
% printYesNo = 1;
% if printYesNo == 1
%     saveFigurePath = '/Users/kevin/SkyDrive/KTH Work/LaTeX Reports/Linear and Nonlinear Hyperbolic Systems/Figures/';
%     set(figure(1), 'PaperPositionMode', 'auto');
%     print('-depsc2', [saveFigurePath ...
%         sprintf('linearizedTime_t_eq_1')]);
% end