%% Coupled Chebfun-FD scheme for time dependent one-dimensional Richards equation

tic
clear all; close all;
T= 0.1;
y1=0; y2=1; yy=0:0.025:1;
D=[0 T y1 y2];
x = chebfun2(@(x,y) x,D);
y = chebfun2(@(x,y) y,D);
E=chebfun2(@(x,y) -x.*y.*(1-y)-0.5,[0 T 0 1]);
phi0 = @(u) diffx(u)./((1-u).^2)-u.^2.*diffy(u,2)-2*u.*diffy(u).*(diffy(u)+1);
phis=phi0(E);                    % source term as chebfun
u0=E(0,:); p0=u0(yy);            % exact solution given by the chebfun2 'E' evaluated at 't=0' and 'y'
t=0; iT=1; nrmv = 1; rezu=1;     % Initialize variables
L=1e3; maxr= 0.5; S=1e6;
corr=zeros(1,S); err_u=zeros(1,S); rez_u=zeros(1,S);
theta = @(p0)  1./(1-p0); tht0 = theta(p0);
pa=p0;
while t<=T
    k=1; s=1; nrmv = 1;
    while nrmv > 1e-5
        [p,tht,dt] = R2D(p0,tht0,t,phis,L,maxr);
        p0=p; tht0=tht;
        nt=41;
        st = linspace(0,1,nt)';
        u=chebfun.spline(st,p);              % convert 1-dim vector to chebfun
        v=u-u0; u0=u;
        nrmv = norm(v);
        if mod(k,10)==0 && iT==1
            u_ref=E(t,:);
            q(:,iT)=p;
            u2=chebfun2(q,[0 T 0 1],'equi'); % convert 2-dim array to chebfun2
            f=phi0(u2)-phis;                 % modified equation with manufactured solution
            corr(s)=nrmv;                    % correction norm
            err_u(s)=norm(u-u_ref);          % error norm
            rez_u(s)=norm(f);                % residual norm
            fprintf('s = %d  ||v||= %2.2e  ||phi(u)||= %2.2e  ||u-u_ref||= %2.2e  \n',s,nrmv,rezu,norm(u-u_ref));
            s=s+1;
        end
        k=k+1; pa=p;
    end
    t=t+dt; iT=iT+1;
    if mod(iT,1)==0
        fprintf('t = %d  iT = %d \n',t,iT);
        di=2;
        xs=1:s-1;
        figure(1);  hold on
        semilogy(xs(1:di:s-1), corr(1:di:s-1),'m*-','LineWidth',1);
        semilogy(xs(1:di:s-1), err_u(1:di:s-1),'gx-','LineWidth',1);
        semilogy(xs(1:di:s-1), rez_u(1:di:s-1),'b+-','LineWidth',1);
        xlabel('Iteration count $s$','Interpreter','latex');
        legend('Corrections = $\|v\|$','Errors = $\|u-u^*\|$','Residuals = $\|\Phi(u)\|$','Interpreter','Latex');
        legend('Box','off','Location','best');
        set(gca,'yscale','log'); box on;
    end
end
figure(2); dk=1; uu=u(yy);
plot(uu);
zlabel('$\psi(z)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
toc
%% FD function ---> solution at fixed time level
function [p,tht,dt] = R2D(p0,tht0,t,phis,L,maxr)
%% Grid Initialization
I = 41;
a=0; b=1;
dx = (b-a)/(I-1);
x = a:dx:b;
%% Parameterization functions
theta = @(c)  1./(1-c);
K = @(c) c.^2; %
%% Initial Conditions
c0=p0; c = c0;
%% solution
tht = theta(c);
DK=K(c);
D=(DK(1:I-1)+DK(2:I))/2;
dt=L*dx^2*maxr/max(D)/2;
F=phis(t,:); fs=F(x);
r=D*maxr/(2*max(D));
rloc=1-(r(1:I-2)+r(2:I-1));
cc(2:I-1)=rloc.*c(2:I-1)+r(2:I-1).*c(3:I)+r(1:I-2).*c(1:I-2);
%% boundary conditions
%%%% BC_Left/Right
cc(1)=p0(1);
cc(I)=p0(I);
%% Source term
dtht=(tht0-tht)/L;
f=diff(r)*dx+dtht(2:I-1)+fs(2:I-1)*dt/L;
cc(2:I-1)=cc(2:I-1)+f;
c=cc;
p=c; tht=theta(c);
end
% Elapsed time is 50.115578 seconds.
