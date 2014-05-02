%% advect - Program to solve the advection equation 
% using the various hyperbolic PDE schemes
clear all;  help advect;  % Clear memory and print header

%% * Select numerical parameters (time step, grid spacing, etc.).

N = 200
L = 1.;     % System size
H = 1;
g = 9.8;
dx = L/N;    % Grid spacing
a = 0.1;
c = g*sqrt(H);      % Wave speed
tau = dx/c;
% coeff = -c*tau/(2.*dx);  % Coefficient used by all schemes
nStep = 5*L/(c*tau);

%% * Set initial and boundary conditions.
sigma = 0.1;              % Width of the Gaussian pulse
k_wave = pi/sigma;        % Wave number of the cosine
x = ((1:N)-1/2)*dx;  % Coordinates of grid points
% Initial condition is a Gaussian-cosine pulse
h = H + a*exp(-(x-L/2).^2/(2*sigma^2)); 
m = c*ones(1,N);
flux = 1/2*g*h.^2 + m.^2/h;
% Use periodic boundary conditions
ip = (2:(N-1))+1;  %ip(N) = N-1;   % ip = i+1 with periodic b.c.
im = (2:(N-1))-1; % im(1) = 2;   % im = i-1 with periodic b.c.

%% * Initialize plotting variables.
iplot = 1;          % Plot counter
hplot(:,1) = h(:);  % Record the initial state
mplot(:,1) = m(:);
tplot(1) = 0;       % Record the initial time (t=0)
nplots = 50;        % Desired number of plots
plotStep = nStep/nplots; % Number of steps between plots

%% * Loop over desired number of steps.
for iStep=1:nStep  %% MAIN LOOP %%
    h(1) = h(2);
    m(1) = -m(2);
    h(N) = h(N-1);
    m(N) = -m(N-1);
    hOld = h;
    mOld = m;
    m(2:(N-1)) = .5*(mOld(ip)+mOld(im)) + -tau/(2*dx)*(flux(ip)-flux(im));
    h(2:(N-1)) = .5*(hOld(ip)+hOld(im)) + -tau/(2*dx)*(mOld(ip)-mOld(im));
    flux = 1/2*g*h.^2 + m.^2/h;
  %* Periodically record a(t) for plotting.
  if( rem(iStep,plotStep) < 1 )  % Every plot_iter steps record 
    iplot = iplot+1;
    hplot(:,iplot) = h(:);       % Record a(i) for ploting
    mplot(:,iplot) = m(:);
    tplot(iplot) = tau*iStep;
    fprintf('%g out of %g steps completed\n',iStep,nStep);
  end
end

%% * Plot the initial and final states.
figure(1); clf;  % Clear figure 1 window and bring forward
plot(x,hplot(:,1),'-',x,h,'--');
legend('Initial  ','Final');
xlabel('x');  ylabel('a(x,t)');
pause(1);    % Pause 1 second between plots

%% * Plot the wave amplitude versus position and time
% figure(2); clf;  % Clear figure 2 window and bring forward
% mesh(tplot,x,hplot);
% ylabel('Position');  xlabel('Time'); zlabel('Amplitude');
% view([-70 50]);  % Better view from this angle

%% Plot
figure(3); clf
for i = 1:nplots
   clf;
   plot(x,hplot(:,i),'-',x,h,'--');
   pause(0.1);
end
%% Plot m
figure(4); clf
for i = 1:nplots
   clf;
   plot(x,mplot(:,i),'-',x,m,'--');
   pause(0.1);
end