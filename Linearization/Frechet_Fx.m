%% Frechet linearization & L-schemes - fixed-point problem phi(u)=0;

tic
clear all; close all; 
cheboppref.setDefaults('plotting','off','display','iter')
% iop=1; % explicit L-scheme
% iop=2; % implicit L-scheme
iop=3; % Frechet linearization
D=[-1 1]; xx=-1:0.02:1; 
x = chebfun('x',D);
A = chebop(D);
a=10; 
u_ref = chebfun(0,D);                              % Exact solution
% u0=chebfun(0.5,D);                               % Initial guess
u0=chebfun('5*sin(pi*x/2)',D);                     % Initial guess
phi = @(u) u+u.^2+a*u.^3;                          % operator of the nonlinear problem
dF = @(u) 3*a*u.^2+2*u+1;
nrmv = 1; y = 0;                
corr=zeros(1,20); rez_v=zeros(1,20); rez_u=zeros(1,20);
rez_uinit=zeros(1,20); err_u=zeros(1,20);
u=u0; 
s=0;
plot3(x,chebfun(y,D),u0,'LineWidth',0.75), hold on
while nrmv > 1e-6
    if iop==1     % (Elapsed time is 2.96 minutes.)
        L=200; % use L=5; for u0=chebfun(0.5,D) % 
        u=u0-sign(dF(u))*phi(u0)/L;                % explicit L-scheme 
        v = u-u0;
        if mod(s,100)==0; fprintf('s = %d  ||v||= %2.2e  \n',s,norm(v)); end
    elseif iop==2 % (Elapsed time is 30.75 seconds.)
        L=5;
        A.op = @(u) L*(u-u0)+sign(dF(u0)).*(u+u0.*u+a*u0.^2.*u);  % implicit L-scheme 
        u=A\0;  
        v = u-u0;
    elseif iop==3 % (Elapsed time is 9.41 seconds.)
        A.op = @(v) (3*a*u.^2+2*u+1)*v+phi(u);     % Newton-Frechet in function space
        v = A\0;                                   
        u = u+v;
    end
    nrmv = norm(v); y = y+norm(v);
    u0=u;
    plot3(x,chebfun(y,D),u)
    xlabel('$x$','Interpreter','latex');
    ylabel('$y=\sum{\|v_k\|}$','Interpreter','latex');
    zlabel('$u_k$','Interpreter','latex'); box on
    s=s+1;
    corr(s)=nrmv;
    err_u(s)=norm(u-u_ref);
    rez_u(s)=norm(phi(u));
end
figure;
semilogy(corr,'m*-','LineWidth',1); hold on
semilogy(rez_u,'gs-','LineWidth',1);
semilogy(err_u,'b+-','LineWidth',1);
xlabel('Iteration count $k$','Interpreter','latex');
legend('Corrections = $\|v\|$','Residuals = $\|\Phi(u)\|$','Errors = $\|u-u^*\|$','Interpreter','latex');
legend('Box','off','Location','best')
figure; dk=1; uu=u(xx);
plot(xx(1:dk:end),uu(1:dk:end),'-+r');
xlabel('$x$','Interpreter','latex'); ylabel('$u(x)$','Interpreter','latex');

toc

%%
