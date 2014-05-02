%% advect - Program to solve the advection equation 
% using the various hyperbolic PDE schemes
clear all;  help advect;  % Clear memory and print header

%% * Select numerical parameters (time step, grid spacing, etc.).

N = 100
L = 1.;     % System size
H = 1;
dx = L/N;    % Grid spacing
c = 1;      % Wave speed
tau = dx/c;
coeff = -c*tau/(2.*dx);  % Coefficient used by all schemes
nStep = L/(c*tau);

%% * Set initial and boundary conditions.
sigma = 0.1;              % Width of the Gaussian pulse
k_wave = pi/sigma;        % Wave number of the cosine
x = ((1:N)-1/2)*dx - L/2;  % Coordinates of grid points
% Initial condition is a Gaussian-cosine pulse
h = H + exp(-x.^2/(2*sigma^2)); 
% Use periodic boundary conditions
ip(1:(N-1)) = 2:N;  ip(N) = 1;   % ip = i+1 with periodic b.c.
im(2:N) = 1:(N-1);  im(1) = N;   % im = i-1 with periodic b.c.

%% * Initialize plotting variables.
iplot = 1;          % Plot counter
hplot(:,1) = h(:);  % Record the initial state
tplot(1) = 0;       % Record the initial time (t=0)
nplots = 50;        % Desired number of plots
plotStep = nStep/nplots; % Number of steps between plots

%% * Loop over desired number of steps.
for iStep=1:nStep  %% MAIN LOOP %%
    h(1:N) = .5*(h(ip)+h(im)) + coeff*(h(ip)-h(im));
  %* Periodically record a(t) for plotting.
  if( rem(iStep,plotStep) < 1 )  % Every plot_iter steps record 
    iplot = iplot+1;
    hplot(:,iplot) = h(:);       % Record a(i) for ploting
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
figure(2); clf;  % Clear figure 2 window and bring forward
mesh(tplot,x,hplot);
ylabel('Position');  xlabel('Time'); zlabel('Amplitude');
view([-70 50]);  % Better view from this angle

%% Plot
figure(3); clf
for i = 1:nplots
   plot(x,hplot(:,i),'-',x,h,'--');
   pause(0.2);
   clf;
end