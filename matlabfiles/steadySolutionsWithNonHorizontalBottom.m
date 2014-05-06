
%% advect linearization
clear all;
clc;
close all;
addpath(['/Users/kevin/SkyDrive/KTH Work/LaTeX Reports',...
    '/HW4-High Resolution shock-capturing methods/matlabfiles/']);
%% Parameters
N = 80;
L = 10.;
dx = L/N; % Grid spacing
H = 1;
g = 9.8;
c = sqrt(g*H); % Wave speed
tau = .8*dx/c; % Time Step
coeff = -tau/(2*dx);
nStep = L/(c*tau);
finalTime = nStep*tau;
nCells = 1:N;
ghostCellOneSide = 2;
cellsPlusGhost = 1:(N+2*ghostCellOneSide);
%% Boundary Conditions
w = 0.1*L; % Width
x = (nCells-1/2)*dx; % Grid points
a = 1/10*H;
h = H + a*exp(-(x-L/2).^2/(w^2));
h = padarray(h',ghostCellOneSide)';
% Boundary Conditions
% m=h;
m = 0*a*exp(-(x-L/2).^2/(w^2))*c;
m = padarray(m',ghostCellOneSide)';
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
i = ((1+ghostCellOneSide):(N+ghostCellOneSide));
im = i-1;
ip = i+1;
%% Record the initial state
hplot(:,1) = h(:);
mplot(:,1) = m(:);
flux = (1/2*g*h.^2 + (m.^2)./h);
r = L/6;
B0 = H/10;
for lengthOfX = 1:length(x)
    if abs(x(lengthOfX)-L/2)<r
        B(lengthOfX) = B0*cos((pi*(x(lengthOfX)-L/2))/(2*r));
    else
        B(lengthOfX)=0;
    end
end
B = padarray(B',ghostCellOneSide)';
B(1) = B(4);
B(2) = B(3);
B(N+4) = B(N+1);
B(N+3) = B(N+2);
    %% Loop
    for iStep=1:(round(nStep))
        hOld=h;
        mOld=m;
        flux = (1/2*g*hOld.^2 + (mOld.^2)./hOld);
        h(i) = .5*(hOld(i+1)+ hOld(i-1)) - abs(coeff)*(   m(i+1)	-m(i-1)) + tau*B(i);
        m(i) = .5*(mOld(i+1)+ mOld(i-1)) - abs(coeff)*(flux(i+1)-flux(i-1)) + tau*B(i);
        h(1) = h(4);
        h(2) = h(3);
        m(1) = m(4)*(-1);
        m(2) = m(3)*(-1);
        h(N+4) = h(N+1);
        h(N+3) = h(N+2);
        m(N+4) = m(N+1)*(-1);
        m(N+3) = m(N+2)*(-1);
        hplot(:,(iStep+1)) = h(:);
        mplot(:,(iStep+1)) = m(:);
    end
    
    %% Plot
    % figure(1); clf;
    % plot(x,hplot(i,1),'-',x,hplot(i,nStep),'--');
    % legend('t=0','t=1');
    % xlabel('x');  ylabel('h(x,t)');
    % title(sprintf('$h(x,t)$ $\\Delta t=$%0.2g,',...
    %     ' $\\epsilon =$ %g',tau,a),...
    %     'Interpreter','latex')
    %% h plot1
    figure(2)
    for iPlotting = 1:2:nStep
        clf;
        plot(x,hplot(i,1))
        hold on;
        plot(x,hplot(i,iPlotting),'-');
        finalB = plot(x,B(i),'-.','color','r');
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
    %% Plot for 2.1
    figure(3)
    initialH = plot(x,hplot(i,1),'-','color','b');
    hold on;
    finalH = plot(x,hplot(i,round(nStep)),'--');
    finalB = plot(x,B(i),'-.','color','r');
    legend([initialH finalH,finalB],...
        't=0',...
        sprintf('t=%0.2g',finalTime),...
        'B(x)');
    xlabel('x');  ylabel('h(x,t)');
    hold off;
    %%
    printYesNo = 1;
    if printYesNo == 1
        saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work/'...
            'LaTeX Reports/HW4-High Resolution shock-capturing '...
            'methods/Figures/'];
        set(figure(3), 'PaperPositionMode', 'auto');
        print('-depsc2', [saveFigurePath ...
            sprintf('steadySolutionsp%g_n_is_%g_a_%g',H,N,round(a))]);
    end
    
    
