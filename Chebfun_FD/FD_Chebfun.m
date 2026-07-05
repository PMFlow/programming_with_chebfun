%% Chebfun solution for unsteady 1D Richards equation with FD time stepping
clear all; close all;
tic
z1=0; z2=1; D=[z1 z2];                   % domain - steady-state problem
z = chebfun('z',D);       
zz=z1:0.1:z2;
N = chebop(D);
T=1; t=0; dt=0.1; p=zeros(11,10);
Dt=[0 T D];                              % domain - unsteady problem
x = chebfun2(@(x,y) x,Dt);
y = chebfun2(@(x,y) y,Dt);
Et=chebfun2(@(x,y) -x.*y.*(1-y)-0.5,Dt); % exact solution
% Richards equation: phi0=phis0
phi0 = @(u) diffx(u)./((1-u).^2)-u.^2.*diffy(u,2)-2*u.*diffy(u).*(diffy(u)-1);
phis0=phi0(Et);
u0=chebfun(-0.5,D);
phi = @(u) u^2*diff(u,2)+2*u*diff(u)*(diff(u)-1);
for i=1:T/dt
    t=t+dt;
    if mod(i,100)==0
        fprintf('i = %d  t = %2.2e  \n',i,t);
    end
    E=chebfun(@(z) -t.*z.*(1-z)-0.5,D);
    phis=phi(E);
    N.op=@(u) -u+u0+(1-u)^2*dt*(u^2*diff(u,2)+2*u*diff(u)*(diff(u)-1)-phis);
    N.lbc=@(u) u-E(z1);
    N.rbc=@(u) u-E(z2);
    u=N\0; u0=u;
    p(:,i)=u(z1:0.1:z2);
end
w=chebfun2(p,[0 1 0 1],'equi'); 
figure; plot(w);
xlabel('$t$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); zlabel('$u(z,t)$','Interpreter','latex'); box on;
figure; plot(w-Et);
xlabel('$t$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); zlabel('$u(z,t)-u^{*}(z,t)$','Interpreter','latex'); box on;
f=phi0(w)-phis0;
rez=norm(f); err=norm(w-Et);
fprintf('Residual = %2.2e  Error = %2.2e  \n',rez,err);
toc
%% dt=0.1; residual = 0.0606; norm(w-Et) = 0.0088; Elapsed time is 7.432385 seconds.
