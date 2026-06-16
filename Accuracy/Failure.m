%%   Failure of modeling saturated-unsaturated transition in Chebfun

close all; clear all;
tic
x1=0; x2=2; D=[x1 x2];
x = chebfun('x',D);
load IC_Richy_sat_101
p=c;
I = length(p);
E=chebfun(p',D,'equi');
c=10;
phi1 = @(u) diff(u,2)+c*diff(u).*(diff(u)+1); % unsaturated 
phi2 = @(u) diff(u,2);                        % saturated
[p,pE,eps] = FDL_scheme(E,I);
u=chebfun(p',D,'equi');
pm=max(0,-p);                                 % negative part of p (unsaturated)
pp=max(0,p);                                  % positive part of p (saturated)
up=chebfun(pp',D,'equi');
um=chebfun(pm',D,'equi');
phi = @(u,up,um) (1-sign(u))/2.*phi1(-um)+(1+sign(u))/2.*phi2(up);
rezu=norm(phi(u,up,um));
disp('Errors');
L2_p = norm(p-pE);
R_p = rezu;
fprintf('FD error ||p-p*|| = %0.2e \n',L2_p)
fprintf('Residual ||Phi(u)|| = %0.2e  \n',R_p)
figure; hold on
xx=0.02:0.02:2;
plot(xx(1:2:end),up(xx(1:2:end)),'+-')
plot(xx(1:2:end),um(xx(1:2:end)),'*-')
plot(xx(1:2:end),u(xx(1:2:end)),'k-','Linewidth',1)
xlabel('$x$','Interpreter','latex'); box on;
legend('$u^+(x)$','$u^-(x)$','$u(x)=u^+(x)-u^-(x)$','Interpreter','latex','box','off','location','best');

figure; hold on
plot(phi(u,up,um),'-r','Linewidth',1)
xlabel('$x$','Interpreter','latex'); box on;
legend('$\Phi^-(-u^-)+\Phi^+(u^+)$','Interpreter','latex','box','off','location','best');

figure; hold on
plot(up,'k','Linewidth',1)
plot(diff(up),'-b','Linewidth',1)
xlabel('$x$','Interpreter','latex'); box on;
legend('$u^+(x)$','$\frac{\partial u^+}{\partial x}$','Interpreter','latex','box','off','location','best');

figure; hold on
plot((1+sign(u))/2.*diff(up,2),'-r','Linewidth',1)
xlabel('$x$','Interpreter','latex'); box on;
legend('$\frac{1+\mbox{sign}(u)}{2}\frac{\partial^2 u^+}{\partial x^2}$','Interpreter','latex','box','off','location','best');
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,pE,eps] = FDL_scheme(E,I)
%%   Grid Initialization
z1 = 0 ; z2 = 2;
dz = (z2-z1)/(I-1);
z = z1:dz:z2;
S = 1e7;
Tolerance=1e-6;
maxr=1;
L=1; % Stabilization parameter
%%  Coefficients
Ksat = 2.77*10^-6; a=Ksat;
alpha=10; c=alpha;
K = @(p) a*exp(+c*p);
%% Reference solution and boundary conditions
pE = E(z); pT=pE(1); pB=pE(I);
p=pE; 
%% Solution
pp=zeros(1,I); eps=zeros(1,S); 
dt=L*dz^2*maxr/2;
for s=1:S
    %% FD scheme:
    DK=K(p);
    D=(DK(1:I-1)+DK(2:I))/2;
    r=dt*D/dz^2/L;
    rloc=1-(r(1:I-2)+r(2:I-1));
    pp(2:I-1)=rloc.*p(2:I-1)+r(2:I-1).*p(3:I)+r(1:I-2).*p(1:I-2);
    pp(1)=pT; pp(I)=pB;
    %% Source term (from L-scheme)
    pp(2:I-1)=pp(2:I-1)+diff(r)*dz;
    p=pp;
    %% Convergence criterion
    eps(s)=norm(p-pE);
    if eps(s)<= Tolerance
        break
    end
end
end