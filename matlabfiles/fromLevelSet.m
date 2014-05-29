%% Level Set Method
% 
%%
function [dt,dx,phiplot,X,Y,tplot,figure1] = levelSetMethod(dtfraction,myTextLabel)
global n xc yc r;
% Long description
if nargin < 2
    n=100;
    dtfraction = 0.1;
    myTextLabel = '0.1';
end
%% Constants
xc = 0;
yc = -.6;
r = 0.3;
L = 2;
dx = L/(n-1);
dy = dx;
% Make grid of x and y values
[X, Y] = meshgrid(-L/2:dx:L/2,-L/2:dy:L/2);
tFinal = 1.6;
x=X(1,:);
dt = dtfraction*dx;
tSteps = ceil(tFinal/dt);
ia = 1:n;
im = [n,1:n-1];
ip = [2:n,1];
% Signed Distance Function
w=.2;
phi = exp(-x.^2/(w^2));
plot(x,phi)
%% Define velocity field
u=-cos(pi*(X+0.5)).*sin(3*pi/8*Y);
v=sin(pi*(X+0.5)).*cos(3*pi/8*Y);
up = u;
um = u;
vp = v;
vm = v;
up(up<=0) = 0;
um(um>=0) = 0;
vp(vp<=0) = 0;
vm(vm>=0) = 0;

%% Plotting variables
nplots = 20;        % Desired number of plots
plotStep = round(tSteps/nplots); % Number of steps between plots
iplot = 1;
tplot(1) = 0;
phinew = phi;
phiplot(:,:,1) = phi;
fx = zeros(n);
fy = zeros(n);
% tplot = zeros(1,nplots);
myHandle = waitbar(0,'Initializing waitbar...');
tic;
%% Main loop
for tn = 1:tSteps
    for i = 2:n-1
            wxm = phi(i) - phi(i-1);	% x backward difference
            wxp = phi(i+1) - phi(i); 	% x forward difference
            fx(i) = vp(i)*wxm + vm(i)*wxp;
            fy(i) = up(i)*wym + um(i)*wyp;
            phinew(i) = phi(i)-dt*( fx(i)/dx + fy(i)/dx );
    end
    % Set Boundary Conditions (no flux)
    phinew(1,:) = phinew(2,:);
    phinew(:,1) = phinew(:,2);
    phinew(n,:) = phinew(n-1,:);
    phinew(:,n) = phinew(:,n-1);
    phi = phinew;
    phiplot(:,:,tn+1) = phi; 
    tplot(tn) = dt*tn;
        %% Here's the progress bar code
    time=toc;
    Perc=tn/tSteps;
    Trem=time/Perc-time; %Calculate the time remaining
    Hrs=floor(Trem/3600);Min=floor((Trem-Hrs*3600)/60);
    waitbar(Perc,myHandle,[sprintf('%0.1f',Perc*100) '%, '...
        sprintf('%03.0f',Hrs) ':'...
        sprintf('%02.0f',Min) ':'...
        sprintf('%02.0f',rem(Trem,60)) ' remaining Level Set Method']);
end
%% Plot for print


qp = 1:round(n/20):n;
totTime=length(tplot);
n1 = 1;
n2 = round(1/4*totTime);
n3 = round(.5*totTime);
n4 = round(3/4*totTime);
n5 = totTime;

t1 = tplot(n1);
t2 = tplot(n2);
t3 = tplot(n3);
t4 = tplot(n4);
t5 = tplot(n5);
figure1=figure('Units', 'pixels', ...
    'Position', [100 100 500 375]);hold on;
y1_plot = contour(X,Y,phiplot(:,:,n1),[0,0],'r'); % initial state
y2_plot = contour(X,Y,phiplot(:,:,n2),[0,0],'m'); % second state
y3_plot = contour(X,Y,phiplot(:,:,n3),[0,0],'color',[.2 .5 .9]); % third state
y4_plot = contour(X,Y,phiplot(:,:,n4),[0,0],'b'); % fourth state
y5_plot = contour(X,Y,phiplot(:,:,n5),[0,0],'color',[.3 .3 .3]); % final state
quivP = quiver(X(qp,qp),Y(qp,qp),u(qp,qp),v(qp,qp),...
    'color',[.5 .5 .5]);
axis([-L/2 L/2 -L/2 L/2])
axis('square')
hXLabel = xlabel('x');
hYLabel = ylabel('y');
% Create axes
hText = text(-0.95,0,...
    [sprintf('N = %g\n\\Deltat = ',n) myTextLabel sprintf(' = %4.4f\n\\Deltax = L/(N-1) = %4.4f',dt,dx)],...
    'EdgeColor',[0 0 0],...
    'BackgroundColor',[1 1 1]);
hLegend = legend(gca,...
    sprintf('t = %0.2f s',t1),...
    sprintf('t = %0.2f s',t2),...
    sprintf('t = %0.2f s',t3),...
    sprintf('t = %0.2f s',t4),...
    sprintf('t = %0.2f s',t5),'location','SouthEast');
hTitle = title(['Contour Plot of '...
    sprintf('$\\mathbf{\\phi}$ with velocity vectors')],...
    'Interpreter','latex');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel, hText], ...
    'FontName'   , 'AvantGarde');
set([hLegend, gca]             , ...
    'FontSize'   , 8           );
set([hXLabel, hYLabel, hText]  , ...
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
    'YGrid'     ,'on',...
    'XColor'      , [.3 .3 .3]    , ...
    'YColor'      , [.3 .3 .3]    , ...
    'LineWidth'   , 1             );
hold off;
%%
% savePath = ['/Users/kevin/SkyDrive/'...
%     'KTH Work/Period 3 2014/DN2255/Homework/5/matlab/mat/'];
% save([savePath sprintf('dt_1-0050e-4',dt)],...
%     'X','Y','u','v','x','y','tSteps','phi','dt')
% save([savePath sprintf('dt1-0050e-4phiplot',dt)],...
%     'X','Y','u','v','x','y','tSteps','phiplot','dt')
end % function
