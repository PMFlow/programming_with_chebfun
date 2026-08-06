%% L-scheme/Frechet switching solutions with Basha model for inceeasing uptake depth h  (Eqs. 16 and 18)

close all; clear all;
cheboppref.setDefaults('factory') % default
tic
n=3;                                    % exponent in the constitutive law K(\psi)=1/(1+\psi^n)
q=0.1;                                  % constant flux on top boundary
b=1;                                    % constant head \psi=b on bottom boundary
x1=0; x2=2; D=[x1 x2];
A = chebop(D);                          % Operator on the interval [0,5]
x = chebfun('x',D);                     % The chebfun "x" on the problem interval
ValueArray = {'o','+','x','*','s'}';
hh=[0 0.5 1 1.5 2];
L=1e2;
for i=1:4
    h=hh(i);
    [F]=step_h(x,h,q);                        % piecewise constant water uptake
    leg{i}=['h =',num2str(h)];
    phi = @(u) (1+u^n).*diff(u,2)-n*u.^(n-1).*diff(u).*(diff(u)+1)+(1+u.^n).^2.*F;
    u0=chebfun(1,D); u = u0;                  % Initial guess
    nrmv = 1;
    s=1;
    while nrmv > 1e-3                                                                              % L-scheme iterations
        A.op = @(u) L*(u-u0)-((1+u0^n)*diff(u,2)-n*(u0^(n-1))*diff(u0)*(diff(u)+1)+(1+u0^n)^2*F);  % Stabilized linear operator
        A.lbc = @(u) diff(u)-(1+u0^n)*q+1;                                                         % Neumann condition at x=x1
        A.rbc = @(u) u-b;                                                                          % Dirichlect condition at x=x2
        u = A\0;                                                                                   % Solve the linearized BVP
        v = u-u0; nrmv = norm(v);                                                                  % Corrections norm
        u0=u;
        s=s+1;
    end
    figure(100)
    plot(u,x,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5); hold on
    set(gca,'Ydir','reverse'); set(gca,'XAxisLocation','top')
    xlabel('$\psi(z)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); ylim([0 2]);
    legend(leg,'Box','off','Location','best')
    u_init=u;
    Function_Frechet_Basha_Fig2(i,u_init);
end
figure(101); legend(leg,'Box','off','Location','best')
toc
function [F]=step_h(x,h,q1)
F_neg = h*q1;
F_pos = 0;
neg = x <= h;
F = F_neg.*neg + F_pos.*~neg;
end

%% function for L-scheme/Frechet switching

function Function_Frechet_Basha_Fig2(i,u_init)
cheboppref.setDefaults('factory')
n=3;                             % exponent in the constitutive law K(\psi)=1/(1+\psi^n)
q=0.1;                           % constant flux on top boundary
b=1;                             % constant head \psi=b on bottom boundary
x1=0; x2=2; D=[x1 x2];
A = chebop(D);                   % Operator on the interval [0,5]
x = chebfun('x',D);              % The chebfun "x" on the interval [0,5]
ValueArray = {'o','+','x','*','s'}';
hh=[0 0.5 1 1.5];
h=hh(i);
leg{i}=['h =',num2str(h)];
fprintf('h= %2.2f i= %d ***********************************************\n',h,i);
[F]=step_h(x,h,q);                                                              % piecewise constant water uptake
phi = @(u) (1+u.^3).*diff(u,2)-3*u.^2.*diff(u).*(diff(u)+1)+(1+u.^3).^2.*F;     % ODE part of the nonlinear BVP - n=3
u = u_init;                                                                     % Initial guess
nrmv = 1;
s=0;
S=20; corr=zeros(1,S); rez_v=zeros(1,S); rez_u=zeros(1,S);
while nrmv > 1e-11
    Du = diff(u);
    A.op = @(v) (1+u^3)*diff(v,2)-(6*u^2*diff(u)+3*u^2)*diff(v)+(3*u^2*diff(u,2)-6*u*diff(u)*(diff(u)+1)+3*u^2*F)*v; % Frechet derivative at u - n=3
    A.lbc = @(v) diff(v)-3*u(x1)^2*v*q+Du(x1)-(1+u(x1)^3)*q+1;                                                       % Neumann condition at x=x1 - n=3
    A.rbc = @(v) v+u(x2)-1;                                                                                          % Dirichlect condition at x=x2
    A.rbc = @(v) v+u(x2)-1;         % Dirichlect condition at x=x2
    v = A\(-phi(u));                % Solve the linearized BVP
    nrmv = norm(v);                 % 2-norm of Newton update
    u = u+v;
    s=s+1;
    corr(s)=nrmv;
    rez_v(s)=norm(A.op(v));
    rez=norm(phi(u));
    rez_u(s)=rez;
end
figure(i);
semilogy(corr,'m*-','LineWidth',1); hold on
semilogy(rez_u,'b+-','LineWidth',1);
semilogy(rez_v,'gx-','LineWidth',1);
title(leg(i));
xlabel('Iteration count $s$','Interpreter','latex');
legend('Corrections = $||v||$','Residuals u = $||\Phi(u)||$','Residuals v = $||A.op(v)||$','Interpreter','latex');
legend('Box','off','Location','best')
figure(101)
plot(u,x,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5); hold on
ylabel('$z$','Interpreter','latex'),xlabel('$\psi(z)$','Interpreter','latex'); ylim([0 2]);
set(gca,'Ydir','reverse'); set(gca,'XAxisLocation','top')
fprintf('op_residual = %2.2e \n',rez_u(s));
end
%% Elapsed time is 53.495374 seconds.