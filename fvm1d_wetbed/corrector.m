function [uhnew, hnew, unew]=corrector(nc,hold,uhold,F,S,dx,dt)

for i=1:nc,
    hnew(i)=hold(i)-dt/dx*(F(i+1,1)-F(i,1));
    uhnew(i)=uhold(i)-dt/dx*(F(i+1,2)-F(i,2))-dt*S(i);
    unew(i)=uhnew(i)/hnew(i); %update for dry bed cases
    %Add vfr for dry bed cases to compute eta from h
end
