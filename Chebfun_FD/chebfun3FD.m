%% %% Coupled Chebfun-FD scheme for time dependent two-dimensional Richards equation

tic
clear all; close all;
T=0.03;
x1=0; x2=1; xx=0:0.0125:1;
y1=0; y2=1; yy=0:0.0125:1;
[X,Y] = meshgrid(xx,yy);
D=[0 T x1 x2 y1 y2];
x = chebfun3(@(x,y,z) x,D);
y = chebfun3(@(x,y,z) y,D);
z = chebfun3(@(x,y,z) z,D);
E=chebfun3(@(x,y,z) -x.*y.*(1-y).*z.*(1-z)-0.5,[0 T 0 1 0 1]);
phi0 = @(r) diffx(r)./((1-r).^2)-r.^2.*(diffy(r,2)+diffz(r,2))-2*r.*(diffy(r).*diffy(r)+diffz(r).*(diffz(r)+1));
phis=phi0(E);                            % source term as chebfun
u0=E(0,:,:); p0=u0(X,Y);                 % exact solution given by the chebfun3 'E' evaluated at 't=0' and '(X,Y)'
t=0; iT=1; nrmv = 1; rezu=1;             % Initialize variables
L=3e3; maxr= 1; S=1e6;
corr=zeros(1,S); err_u=zeros(1,S); rez_u=zeros(1,S);
theta = @(p0)  1./(1-p0); tht0 = theta(p0);
pa=p0;
while t<=T
    k=1; s=1; nrmv = 1;
    while nrmv > 1e-7
        [p,tht,dt] = R3D(p0,tht0,t,phis,L,maxr);
        p0=p; tht0=tht;
        u=chebfun2(p,[0 1 0 1],'equi');              % convert 2-dim array to chebfun2
        v=u-u0; u0=u;
        nrmv = norm(v);
        if mod(k,50)==0 && round(T-t)==0
            q(iT,:,:)=p;
            u_ref=E(t,:,:);
            u3=chebfun3(q,[0 T 0 1 0 1],'equi');     % convert 3-dim array to chebfun3
            f=phi0(u3)-phis;                         % modified equation with manufactured solution
            corr(s)=nrmv;                            % correction norm
            err_u(s)=norm(u-u_ref);                  % error norm
            rez_u(s)=norm(f);                        % residual norm
            fprintf('s = %d  ||v||= %2.2e  ||phi(u)||= %2.2e  ||u-u_ref||= %2.2e  \n',s,nrmv,rezu,norm(u-u_ref));
            s=s+1;
        end
        k=k+1; pa=p;
    end
    t=t+dt; iT=iT+1;
    if mod(iT,1)==0
        fprintf('t = %d  iT = %d \n',t,iT);
        di=10;
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
figure(2); dk=1; uu=u(X,Y);
surf(uu);
zlabel('$\psi(z)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');xlabel('$x$','Interpreter','latex');
toc
%% FD function ---> solution at fixed time level
function [p,tht,dt] = R3D(p0,tht0,t,phis,L,maxr)
%%   Grid Initialization
I = 81; J = I;
a=0; b=1;
c=0; d=1;
dx = (b-a)/(I-1);
x = a:dx:b;
dy=(d-c)/(J-1);
y=c:dy:d;
[X,Y] = meshgrid(x,y);
%% Parameterization functions
theta = @(c)  1./(1-c);
K = @(c) c.^2; %
%% Initial Conditions
c0=p0; c = c0;
%% solution
cc=zeros(J,I);
tht = theta(c);
DK=K(c);
Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
dt=12*(max(max(Dx))/dx^2/L+max(max(Dy))/dy^2/L); dt=maxr*1/dt;
F=phis(t,:,:); fs=F(X,Y);
rx=dt*Dx/dx^2/L; ry=dt*Dy/dy^2/L;
rloc=1-(rx(:,1:I-2)+rx(:,2:I-1)+ry(1:J-2,:)+ry(2:J-1,:));
cc(2:J-1,2:I-1)=rloc.*c(2:J-1,2:I-1) ...
    +rx(:,1:I-2).*c(2:J-1,1:I-2)+rx(:,2:I-1).*c(2:J-1,3:I) ...
    +ry(1:J-2,:).*c(1:J-2,2:I-1) +ry(2:J-1,:).*c(3:J,2:I-1);
%% boundary conditions
%%%% BCXLeft/Right
cc(:,1)=c0(:,1);
cc(:,I)=c0(:,I);
%%%% BCYBottom/Upper
cc(1,:)=c0(1,:);
cc(J,:)=c0(J,:);
%% Source term
dtht=(tht0-tht)/L;
f=(ry(2:J-1,:)-ry(1:J-2,:))*dy + dtht(2:J-1,2:I-1)...
    + fs(2:J-1,2:I-1)*dt/L;
cc(2:J-1,2:I-1)=cc(2:J-1,2:I-1)+f;
c=cc;
p=c; tht=theta(c);
end
% Elapsed time is 1521.921849 seconds = 25.3654 minutes.
