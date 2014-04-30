function [F, amax0] = fluxes(grav,nf,eta,u,z,dz,deta,du)

%I'm using a 2d solver in 1d, so I'm setting sn=0 and cn=1 for all cases
sn=0;
cn=1;
vl=0;
vr=0;

%Left Boundary : model as wall
hr=eta(1)-0.5*deta(1)-z(1);
ur=u(1)-0.5*du(1);
[fdum, amax]=solver(hr,hr,-ur,ur,vl,vr,sn,cn);
F(1,1)=fdum(1);
F(1,2)=fdum(2);

%Right Boundary : model as wall
hl=eta(nf-1)+0.5*deta(nf-1)-z(nf);
ul=u(nf-1)+0.5*du(nf-1);
[fdum, amax]=solver(hl,hl,ul,-ul,vl,vr,sn,cn);
F(nf,1)=fdum(1);
F(nf,2)=fdum(2);

%Other bc options could be used
%q=12.0;
%hr=2;
%F(1,1)=q;
%F(1,2)=q^2/hr+0.5*grav*hr^2;

%F(nf,1)=hl*ul;
%F(nf,2)=ul^2*hl+0.5*grav*hl^2;



amax0=0;
for i=2:nf-1, %Sweep over faces
    %Variable reconstruction
    hl=eta(i-1)+0.5*deta(i-1)-z(i);
    ul=u(i-1)+0.5*du(i-1);
    hr=eta(i)-0.5*deta(i)-z(i);
    ur=u(i)-0.5*du(i);
    %Call solver
    [fdum, amax]=solver(hl,hr,ul,ur,vl,vr,sn,cn);
    F(i,1)=fdum(1);
    F(i,2)=fdum(2);
    amax0=max([amax0 amax]); %Keep track of max wave speed to check CFL
end








