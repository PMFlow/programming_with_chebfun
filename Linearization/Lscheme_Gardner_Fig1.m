%% L-scheme for homogeneous BVP with Grdner model (Eqs. 16 and 17)

tic
clear all; close all; % clc;
% cheboppref.setDefaults('plotting','off','display','iter','maxIter',500) 
cheboppref.setDefaults('factory'); % default
x1=0; x2=2; D=[x1 x2]; xx=0:0.02:2; L=100;          
b=1;      % constant head \psi=b on bottom
q=0.1;    % constant flux on top
a=1; c=2; % parameters for K(\psi)=a*exp(-c*\psi) 
A = chebop(D); 
x = chebfun('x',D); 
ValueArray = {'o','+','x','*'}';
L=2; 
%% L-scheme: L(u-ua)=OpDiff; L(u-ua)-OpDiff=0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:4
c=0.5*i;
fprintf('c= %2.2f  i= %d ***************************************\n',c,i);
leg{i}=['c=',num2str(c)];
phi = @(u) diff(u,2)-c*diff(u).*(diff(u)+1); % ODE part of the nonlinear BVP
u = x; u0=u;                                 % Initial guess
nrmv = 1; y = 0;                             % Initialize variables
corr=zeros(1,20); rez_v=zeros(1,20); rez_u=zeros(1,20); rez_uinit=zeros(1,20);
s=0; 
figure(i)
plot3(x,chebfun(y,D),u,'LineWidth',0.75), hold on
while nrmv > 1e-11                                             % L-scheme iterations   
    A.op  = @(u) L*(u-u0)-(diff(u,2)-c*diff(u0)*(diff(u)+1));  % Stabilized linear operator - implicit L-scheme
    A.lbc = @(u) a*(diff(u)+1)-q*exp(c*u0);                    % Neumann condition at x=x1
    A.rbc = @(u) u-b;                                          % Dirichlect condition at x=x2
    u = A\0;                                                   % Solve the linearized BVP
    v = u-u0; nrmv = norm(v); y = y+norm(v);                   % 2-norm of Newton update
figure(i)
plot3(x,chebfun(y,D),u)
title(leg(i));
xlabel('$x$','Interpreter','latex'); 
ylabel('$y=\sum{\|v_s\|}$','Interpreter','latex'); 
zlabel('$u_s$','Interpreter','latex'); box on
s=s+1; 
corr(s)=nrmv;
rez_u(s)=norm(phi(u));  
u0=u;
end
figure(4+i); 
semilogy(corr,'m*-','LineWidth',1); hold on
semilogy(rez_u,'b+-','LineWidth',1);
% title(leg(i));
xlabel('Iteration count $k$','Interpreter','latex');
legend('Corrections = $||v||$','Residuals = $||\Phi(u)||$','Interpreter','latex');
legend('Box','off','Location','best')
figure(112)
plot(u,x,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',7); hold on
end
legend(leg,'box','off','Location','best')
set(gca,'yscale','linear');set(gca,'Ydir','reverse');set(gca,'XAxisLocation','top');
xlabel('$\psi(z)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); ylim([0 2]);
toc
% Elapsed time is 24.712937 seconds. 
