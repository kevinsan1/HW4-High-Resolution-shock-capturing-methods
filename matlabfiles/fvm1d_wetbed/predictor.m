function [etap, up]=predictor(nc,eta,h,u,deta,du,dz,zc,dt,dx);

global grav

dh=deta-dz;
for i=1:nc,
    etap(i)=eta(i)-0.5*dt/dx*(h(i)*du(i)+u(i)*dh(i));
    up(i)=u(i)-0.5*dt/dx*(u(i)*du(i) + grav*deta(i));
end