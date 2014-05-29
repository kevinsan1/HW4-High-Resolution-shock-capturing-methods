%FVM1D - FOR WET BED PROBLEMS ONLY
%by Brett F. Sanders (and pieces of code from Scott F. Bradford)
%
%This is a very simple 1D solver of the shallow-water equations
%that uses the Hancock predictor-corrector time-stepping scheme,
%the MUSCL method of slope limiting and variable reconstruction
%and Roe's approximate Riemann solver to compute fluxes.
%
%The code is kept as simple as possible to emphasize the basic
%flow of logic. To account for problems involving a dry bed,
%one needs to add a number of "if" statements to avoid
%division by zero. This makes the code pretty messy and
%therefore these lines have been omitted.
%
%Note that the code can either be run in a first order or
%second order accurate mode. The user can select from 
%several limiters to see how these impact the solution.
%
%To run this program, copy all the .m files into a directory
%and run "fvm1d.m" by either typing "fvm1d" at the matlab
%command prompt or pushing the execute button in the matlab
%text editor.

clear
close all
format long

global grav

%Set up grid
nc=100; %number of cells
nf=nc+1; %number of edges
L=1000; %length of channel
dx=L/nc; %length of cell
x=0:dx:L; %array of edge coordinates
xc=dx/2:dx:L-dx/2; %array of cell center coordinates

%Set up time marching and output interval
dt=0.5; %time step (s)
nt=100; %number of time steps
ntplot=10; %plot interval (number of time steps)

%Define bed elevation at faces z=f(x)
z=zeros(size(x)); %flat bed - can enter own function here, z=f(x)

%Compute bed slope
for i=1:nc,
    dz(i)=z(i+1)-z(i); %dimensions of length
    zc(i)=0.5d0*(z(i+1)+z(i)); %elevation of cell center
end

%Set parameter values
grav=9.806;

%Set attributes of solver
iorder=2; %1=first order scheme, 2=second order scheme
beta=2; %controls limiter used by model
%Notes on limiters
%beta=1 => Minmod
%beta=2 => Superbee
%beta=3 => Fromm scheme, predicts oscillations at sharp fronts
%beta=4 => Van Leer
%beta=5 => Van Albada
%beta=6 => Double Minmod

%Set up initial condition
xo=L/2;
etalo=10;
etaro=1;
ulo=0; 
uro=0;
for i=1:nc,
    if (xc(i) < xo),
        eta(i)=etalo;
        h(i)=eta(i)-zc(i);
        u(i)=ulo;
    else
        eta(i)=etaro;
        h(i)=eta(i)-zc(i);
        u(i)=uro;
    end
end

%Initialize arrays
uh=h.*u;

deta=zeros(size(eta));
du=zeros(size(u));

t=0; %start time
for n=1:nt, %Begin time-marching loop
    if (iorder == 2), %for second order accuracy only
        deta = limiter(nc,beta,eta);
        du = limiter(nc,beta,u);
        [etap, up]=predictor(nc,eta,h,u,deta,du,dz,zc,dt,dx);
        hp=etap-zc; %Update for dry bed cases
        S = grav*hp.*dz/dx; %Source term treatment
        [F, amax] = fluxes(grav,nf,etap,up,z,dz,deta,du); %Compute fluxes
    else
        [F, amax] = fluxes(grav,nf,eta,u,z,dz,deta,du);
        S = grav*h.*dz/dx; %Source term treatment
    end
    [uh, h, u] = corrector(nc,h,uh,F,S,dx,dt);
    eta=h+zc; %compute new free surface height
    e=eta+0.5*u.^2/grav; %compute energy in units of length (head)
    t=t+dt;
    cr=amax*dt/dx; %compute Courant number
    fprintf(1,'%g %d\n',n,cr)
    if (cr > 1) %Stops program if Courant number exceeds one.
        break
    end
    if (mod(n,ntplot) == 0),
        subplot(3,1,1)
        plot(xc,e,'r-',xc,h+zc,'b-',xc,zc,'k-')
        axis([0 L 0 16])
        legend('Energy','Free Surface','Bed')
        subplot(3,1,2)
        plot(xc,u,'b-')
        axis([0 L -1 10]) 
        legend('Velocity')
        subplot(3,1,3)
        plot(x,F(:,1),'b-')
        axis([0 L -1 50]) 
        legend('Discharge') 
        pause(0.1)
    end
end


