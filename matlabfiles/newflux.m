%% advect linearization
clear all;
clc;
close all;
addpath(['/Users/kevin/SkyDrive/KTH Work/Period 3 2014'...
    '/DN2255/Homework/4/HW4-High Resolution shock-capturing'...
    ' methods/matlabfiles/']);
% Parameters
N = 100;
L = 1.;
dx = L/N; % Grid spacing
H = 1;
g = 9.8;
c = sqrt(g*H); % Wave speed
tau = .8*dx/c; % Time Step
coeff = -tau/(2*dx);
nStep = L/(c*tau); % waves
finalTime = nStep*tau;
nCells = 1:N;
x = dx/2:dx:L;
% Boundary Conditions
w = 0.1*L; % Width
a = 1/5*H;
h = a*exp(-(x-L/2).^2/(w^2));
m = 0*h;
% Record the initial state
hplot(:,1) = h(:);
mplot(:,1) = m(:);
f = (1/2*g*h.^2 + (m.^2)./h);
% Change of Variables
q1 = h;
q2 = m;
q = [q1; q2];
z1 = h./sqrt(h);
z2 = m./sqrt(h);
Z = [z1; z2];
qZ = [z1.^2; z2.*z1];
dqdz = [...
    2*z1,   zeros(1,N);
    z2,     z1];
fZ = [...
    z1.*z1;
    z2.^2 + 1/2*g*(z1).^4
    ];
dfdz = [z2 			z1; ...
    2*g*z1.^3	2*z2];
% The average between the endpoints
im = [N,1:(N-1)];
ip = [2:N,1];
myi = 1:N;
p(1) = 1; p(2) = 2;
zpLinearIntegral = 1/2*(Z(p(1),im) + Z(p(1),:));
% For the cubic term
zpCubicTerm = .5*(Z(p(1),im)    + Z(p(1),:)   ).* ...
    (.5*(Z(p(1),im).^2  + Z(p(1),:).^2));
Z1 = 1/2*(Z(p(1),im) + Z(p(1),myi)); % average for p = 1
Z2 = 1/2*(Z(p(2),im) + Z(p(2),myi)); % average for p = 2
% or Z1.*htilde;
htilde = 1/2*(h(im) + h);
zpCubicTerm2 = Z(p(1),:).*htilde;
isequal(zpCubicTerm,zpCubicTerm2);
%
Bm = [ ...
    2*Z1, 	zeros(1,length(Z1)); ...
    Z2, 		Z1];
% %%
Cm = [ ...
    Z2, 				Z1; ...
    2*g*Z1.*htilde, 	2*Z2];
Am = [
    zeros(1,length(Z1)),      ones(1,length(Z1));
    -(Z2./Z1).^2 + g*htilde,  2*Z1./Z1
    ];
uhat = Z2./Z1;
u = m./h;
uhat2 = (sqrt(h(im)).*u(im) + sqrt(h))./...
    (sqrt(h(im)) + sqrt(h));
Am3 = Cm*Bm'
Am2 = [zeros(1,N), ones(1,N);
        -uhat.^2 + g*htilde, 2*uhat];
    isequal(Am2,Am)
chat = sqrt(g*htilde);
lambdahat(1,:) = uhat-chat;
lambdahat(2,:) = uhat+chat;
rhat1 = lambdahat(1,:);
rhat2 = lambdahat(2,:);
Lhat = [(uhat+chat)./(2*chat), -1*ones(1,N)./(2*chat);
    -(uhat-chat)./(2*chat), ones(1,N)./(2*chat)];
D = q(:,myi) - q(:,im);
alpha(1,:) = ( (uhat + chat).*D(1,:) - D(2,:))./(2*chat);
alpha(2,:) = (-(uhat - chat).*D(1,:) - D(2,:))./(2*chat);
W1 = alpha(1,:).*rhat1;
W2 = alpha(2,:).*rhat2;
itest = 2;
fQ1 = q(2,:);
fQ2 = q(2,:).^2./q(1,:) + 1/2*g*(q(1,:)).^2;
htilde = 1/2*(h(im) + h);
chat = sqrt(g*htilde);
for itest = 2:(N-1)
    u(itest) = m(itest)./h(itest);
    uhat2(itest) = (sqrt(h(itest-1)).*u(itest-1) + ...
        sqrt(h(itest)))./...
        (sqrt(h(itest)) + sqrt(h(itest)));
    lambdahat1(itest) = uhat2(itest)-chat(itest);
    lambdahat2(itest) = uhat2(itest)+chat(itest);
    Flux1(itest,:) = 1/2 * (fQ1(itest) + fQ1(itest+1)) - ...
    1/2*(lambdahat(1,itest)*W1 + lambdahat(2,itest)*W2);
    Flux2(itest,:) = 1/2 * (fQ2(itest) + fQ2(itest+1)) - ...
    1/2*(lambdahat(1,itest)*W1 + lambdahat(2,itest)*W2);
end
Flux1(1,:) = Flux1(2,:);
Flux1(N,:) = Flux1(N-1,:);
Flux2(1,:) = Flux2(2,:);
Flux2(N,:) = Flux2(N-1,:);
%%
    Flux = 1/2 * (fQ(itest) + fQ(itest+1)) - ...
    1/2*(lambdahat(1,itest)*W1 + lambdahat(2,itest)*W2);

%%
Flux = 1/2 * (fQ(im) + fQ(ip)) - ...
    1/2*(lambdahat(1,:).*W1 + lambdahat(2,:).*W2);
Q = [q1;q2];
Qn(1:2,1:N,1) = Q;
%%

for ttest = 1:10
   Qn(1:2,2:(N-1),ttest+1) = Q(1:2,2:(N-1)) + ...
       tau/dx * (Flux(3:N) - Flux(1:(N-2)))
   Qn(1:2,1,ttest+1) = Qn(1:2,2,ttest+1);
   Qn(1:2,N,ttest+1) = Qn(1:2,N-1,ttest+1);
   Q = Qn(1:2,1:N,ttest+1);
end

%% Loop
for tstep= 1:(round(nStep))
    hp=h;
    mp=m;
    f = (1/2*g*hp.^2 + (mp.^2)./hp);
    
    for i = 2:(N-1)
        F = 1/2*(f(i) + f(i+1)) - 1/2
        hp(i) = .5*(hp(i+1)+ hp(i-1))...
            - abs(coeff)*(m(i+1) - m(i-1));
        mp(i) = .5*(mp(i+1)+ mp(i-1))...
            - abs(coeff)*(f(i+1) - f(i-1));
    end
    hp(1) = hp(2);
    mp(1) = mp(2)*(-1);
    hp(N) = hp(N-1);
    mp(N) = mp(N-1)*(-1);
    hplot(:,(tstep+1)) = hp(:);
    mplot(:,(tstep+1)) = mp(:);
    h=hp;
    m=mp;
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
%     figure(2)
%     for iPlotting = 1:2:nStep
%         clf;
%         plot(x,hplot(i,1) + B(i)')
%         hold on;
%         plot(x,hplot(i,iPlotting)+B(i)','-');
%         finalB = plot(x,B(i),'-.','color','r');
%         pause(.1)
%     end
% m plot
figure(3)
for ip = 1:8:nStep
    plot(x,hplot(:,1))
    hold on;
    plot(x,hplot(:,ip),'-');
    hold off;
    pause(.1)
end
%% Plot for 2.1
% figure(3)
% initialH = plot(x,hplot(:,1) + B(:)','-','color','b');
% hold on;
% hPlusB.t1 = plot(x,hplot(:,round(nStep/2))+ B(:)'...
%     ,'--','color','k');
% hPlusB.t2 = plot(x,hplot(:,round(nStep))+ B(:)'...
%     ,'--','color','m');
% finalB = plot(x,B(:),'-.','color','r');
% legend([initialH,hPlusB.t1,hPlusB.t2,finalB],...
%     't=0, h(x) + B(x)',...
%     sprintf('t=%0.2g, h(x) + B(x)',tau*round(nStep/2)),...
%     sprintf('t=%0.2g, h(x) + B(x)',tau*round(nStep)),...
%     'B(x)','location','best');
% xlabel('x');  ylabel('h(x,t)');
% hold off;
% %%
% printYesNo = 0;
% if printYesNo == 1
%     saveFigurePath = ['/Users/kevin/SkyDrive/KTH Work'...
%         '/Period 3 2014/DN2255/Homework/4/HW4-High'...
%         ' Resolution shock-capturing methods/Figures/'];
%     set(figure(3), 'PaperPositionMode', 'auto');
%     print('-depsc2', [saveFigurePath ...
%         sprintf('steadySolutionsp%g_n_is_%g_a_%g',H,N,round(a))]);
% end
%
%
