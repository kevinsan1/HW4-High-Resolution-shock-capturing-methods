
%% advect linearization
clear all;
clc;
close all;
addpath(['/Users/kevin/SkyDrive/KTH Work/Period 3 2014'...
    '/DN2255/Homework/4/HW4-High Resolution shock-capturing'...
    ' methods/matlabfiles/']);
%% Parameters
N = 80;
L = 10.;
dx = L/N; % Grid spacing
H = 1;
g = 9.8;
c = sqrt(g*H+1); % Wave speed
tau = .85*dx/c; % Time Step
coeff = -tau/(dx);
nStep = L/(c*tau);
finalTime = nStep*tau;
nCells = 1:N;
ghostCellOneSide = 2;
cellsPlusGhost = 1:(N+2*ghostCellOneSide);
gi = (1:N) + 2;
%% Boundary Conditions
w = 0.1*L; % Width
x = (nCells-1/2)*dx; % Grid points
a = 1/10*H;
h = H + a*exp(-(x-L/2).^2/(w^2));
h = padarray(h,[0,2]); % add 2 ghost cells on both sides
% Boundary Conditions
% m=h;
m = 0*a*exp(-(x-L/2).^2/(w^2))*c;
m = padarray(m,[0,2]); % add 2 ghost cells on both sides
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
u = m./h;
Q = [h;m];
uhat(gi) = (sqrt(h(gi-1)).*u(gi-1) + sqrt(h(gi)))/...
    (sqrt(h(gi-1)) + sqrt(h(gi)));
uhat(N+4) = uhat(N+1)*(-1);
uhat(N+3) = uhat(N+2)*(-1);
uhat(1) = uhat(4)*(-1);
uhat(2) = uhat(3)*(-1);
%% Record the initial state
hplot(:,1) = h(:);
mplot(:,1) = m(:);
fs = (1/2*g*h.^2 + (m.^2)./h);
r = L/6;
B0 = H/10;
for lengthOfX = 1:length(x)
    if abs(x(lengthOfX)-L/2)<r
        B(lengthOfX) = B0*cos((pi*(x(lengthOfX)-L/2))/(2*r));
    else
        B(lengthOfX)=0;
    end
end
B = padarray(B,[0,2]); % add 2 ghost cells on both sides
B(1) = B(4);
B(2) = B(3);
B(N+4) = B(N+1);
B(N+3) = B(N+2);
%% Loop
for iStep=1:(round(nStep))
    fs = (1/2*g*h.^2 + (m.^2)./h);
    f = [m; fs];
    for i = 3:(N+2) % 2 ghost cells on each end 1,2 and N+3,N+4
        Am = [0,1;-uhat(i-1)^2 + g*h(i-1), 2*uhat(i-1)];
        Ap = [0,1;-uhat(i+1)^2 + g*h(i+1), 2*uhat(i+1)];
        Fip = 1/2*(f(:,i) + f(:,i+1)) - ...
            1/2*abs(Ap)*(Q(:,i+1) - Q(:,i));
        Fim = 1/2*(f(:,i) + f(:,i-1)) - ...
            1/2*abs(Am)*(Q(:,i-1) - Q(:,i));
        Qn(:,i) = Q(:,i) + tau/dx*(Fip - Fim);
        flux = (1/2*g*h.^2 + (m.^2)./h);
        hn(i) = .5*(h(i+1)+ h(i-1))...
            - 1/2*abs(coeff)*(m(i+1) - m(i-1));
        mn(i) = .5*(m(i+1)+ m(i-1))...
            - 1/2*abs(coeff)*(flux(i+1) - flux(i-1)...
            + g*hn(i).*(B(i+1) - B(i-1)) );
    end
    hn(1) = hn(4);
    hn(2) = hn(3);
    mn(1) = mn(4)*(-1);
    mn(2) = mn(3)*(-1);
    hn(N+4) = hn(N+1);
    hn(N+3) = hn(N+2);
    mn(N+4) = mn(N+1)*(-1);
    mn(N+3) = mn(N+2)*(-1);
    uhat(gi) = (sqrt(h(gi-1)).*u(gi-1) + sqrt(h(gi)))/...
        (sqrt(h(gi-1)) + sqrt(h(gi)));
    uhat(N+4) = uhat(N+1)*(-1);
    uhat(N+3) = uhat(N+2)*(-1);
    uhat(1) = uhat(4)*(-1);
    uhat(2) = uhat(3)*(-1);
    h=hn;
    m=mn;
    u = m./h;
    Qplot(:,:,iStep+1) = Qn;
    Q = Qn;
    hplot(:,(iStep+1)) = hn(:);
    mplot(:,(iStep+1)) = mn(:);
end
%% Plot
for n = 1:95
    figure(1);clf;
    plot(x,hplot(gi,1),x,hplot(gi,n))
    hold on;
    % hTitle, hXLabel, hYLabel
    hTitle = title(sprintf('Title'));
    hXLabel = xlabel('x');
    hYLabel = ylabel('y');
    hLegend = legend('Initial',sprintf('t=%0.2f',n*tau));
    % Configuration
    set( gca                       , ...
        'FontName'   , 'Helvetica' );
    set([hTitle, hXLabel, hYLabel,hLegend], ...
        'FontName'   , 'AvantGarde');
    set( gca             , ...
        'FontSize'   , 8           );
    set([hXLabel, hYLabel,hLegend]  , ...
        'FontSize'   , 10          );
    set( hTitle                    , ...
        'FontSize'   , 12          , ...
        'FontWeight' , 'bold'      );
    set(gca, ...
        'Box'         , 'off'         , ...
        'TickDir'     , 'out'         , ...
        'TickLength'  , [.02 .02]     , ...
        'XMinorTick'  , 'on'          , ...
        'YMinorTick'  , 'on'          , ...
        'XColor'      , [.3 .3 .3]    , ...
        'YColor'      , [.3 .3 .3]    , ...
        'ZColor'      , [.3 .3 .3]    , ...
        'LineWidth'   , 1             );
    hold off;
    pause(0.001);
end
