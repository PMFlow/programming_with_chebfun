%% L-scheme/Frechet switching solutions with Basha model for different (n,q) - [Basha, 1999, Fig1]  (Eqs. 16 and 18)

close all; clear all;
cheboppref.setDefaults('factory') % default
tic
b=1;                        % constant head \psi=b on bottom
x1=0; x2=2; D=[x1 x2]; xx=0:0.02:2;
A = chebop(D);              % Operator on the interval [0,5]
x = chebfun('x',D);         % The chebfun "x" on the problem interval
ValueArray = {'o','+','x','*'}';
L=150;
for i=1:4
    if i== 1
        n=3; q=0.01;
    elseif i==2
        n=3; q=0.1;
    elseif i==3
        n=7; q=0.01;
    else
        n=7; q=0.1;
    end
    leg{i}=['n=',num2str(n), '   q=',num2str(q)];
    phi = @(u) (1+u^n)*diff(u,2)-n*(u^(n-1))*diff(u)*(diff(u)+1);                      % ODE part of the nonlinear BVP
    u0=chebfun(1,D); u = u0;                                                           % Initial guess
    nrmv = 1; 
    s=1;
    while nrmv > 1e-3                                                                  % L-scheme iterations
        A.op  = @(u) L*(u-u0)-((1+u0^n)*diff(u,2)-n*(u0^(n-1))*diff(u0)*(diff(u)+1));  % Stabilized linear operator
        A.lbc = @(u) diff(u)-(1+u0^n)*q+1;                                             % Neumann condition at x=x1
        A.rbc = @(u) u-b;                                                              % Dirichlect condition at x=x2
        u = A\0;                                                                       % Solve the linearized BVP
        v = u-u0; nrmv = norm(v);                                                      % 2-norm of Newton update
        u0=u;
        s=s+1;
    end
    figure(100)
    plot(u,x,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5); hold on
    set(gca,'Ydir','reverse'); set(gca,'XAxisLocation','top')
    xlabel('$\psi(z)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
    legend(leg,'Box','off','Location','best')
    u_init=u;
    Function_Frechet_Basha_Fig1(i,u_init);
end
figure(101); legend(leg,'Box','off','Location','best')
toc

%% function for L-scheme/Frechet switching

function Function_Frechet_Basha_Fig1(i,u_init)
cheboppref.setDefaults('factory')
x1=0; x2=2; D=[x1 x2]; xx=0:0.02:2; L=0.5;
A = chebop(D);           % Operator on the interval [0,5]
x = chebfun('x',D);      % The chebfun "x" on the interval [0,5]
ValueArray = {'o','+','x','*'}';
if i== 1
    n=3; q=0.01;
elseif i==2
    n=3; q=0.1;
elseif i==3
    n=7; q=0.01;
else
    n=7; q=0.1;
end
fprintf('n= %2.2f  q= %2.2f  i= %d ***************************************\n',n,q,i);
leg{i}=['n=',num2str(n), '   q=',num2str(q)];
if n==3
    phi = @(u) (1+u.^3).*diff(u,2)-3*u.^2.*diff(u).*(diff(u)+1); % ODE part of the nonlinear BVP - n=3
else
    phi = @(u) (1+u.^7).*diff(u,2)-7*u.^6.*diff(u).*(diff(u)+1); % ODE part of the nonlinear BVP - n=7
end
u = u_init;                                                      % Initial guess
nrmv = 1;                                                        % Initialize variables
corr=zeros(1,20); rez_v=zeros(1,20); rez_u=zeros(1,20);
s=0;
while nrmv > 1e-11                          % Newton iterations
    Du = diff(u);                           % Needed to compute the BCs
    if n==3
        A.op = @(v) (1+u^3)*diff(v,2)-(6*u^2*diff(u)+3*u^2)*diff(v)+(3*u^2*diff(u,2)-6*u*diff(u)*(diff(u)+1))*v;     % Frechet differential at u - n=3
        A.lbc = @(v) diff(v)-3*u(x1)^2*v*q+Du(x1)-(1+u(x1)^3)*q+1;                                                   % Neumann condition at x=x1 - n=3
    else
        A.op = @(v) (1+u^7)*diff(v,2)-(14*u^6*diff(u)+7*u^6)*diff(v)+(7*u^6*diff(u,2)-42*u^5*diff(u)*(diff(u)+1))*v; % Frechet differential at u - n=7
        A.lbc = @(v) diff(v)-7*u(x1)^6*v*q+Du(x1)-(1+u(x1)^7)*q+1;                                                   % Neumann condition at x=x1 - n=7
    end
    A.rbc = @(v) v+u(x2)-1;                                                                                          % Dirichlect condition at x=x2
    v = A\(-phi(u));                                                                                                 % Solve the linearized BVP
    nrmv = norm(v);                                                                                                  % 2-norm of Newton update
    u = u+v;
    s=s+1;
    corr(s)=nrmv;
    rez_v(s)=norm(A.op(v));
    rez_u(s)=norm(phi(u));
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
%% Elapsed time is 23.526200 seconds.