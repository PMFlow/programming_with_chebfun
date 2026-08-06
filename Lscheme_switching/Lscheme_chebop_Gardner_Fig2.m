%% L-scheme/chebop switching for inhomogeneous BVP with Gardner model (Eqs. 16 and 17)

tic
clear all; close all;
cheboppref.setDefaults('factory'); % default 
x1=0; x2=2; D=[x1 x2]; xx=0:0.02:2; L=0.5;
b=1;        % constant head \psi=b on bottom
q=0.1;      % constant flux on top
a=1; c=0.5; % parameters for K(\psi)=a*exp(-c*\psi)
A = chebop(D);
x = chebfun('x',D);
ValueArray = {'o','+','x','*','s'}';
hh=[0 0.5 1 1.5];

for i=1:4
    h=hh(i);
    leg{i}=['h =',num2str(h)];
    [F]=step_h(x,h,q);                                        % piecewise constant water uptake
    phi = @(u) diff(u,2)-c*diff(u).*(diff(u)+1)+exp(c*u)*F/a; % ODE part of the nonlinear BVP
    u = x; u0=u;                                              % Initial guess
    nrmv = 1;                                                 
    corr=zeros(1,20); rez_v=zeros(1,20); rez_u=zeros(1,20); rez_uinit=zeros(1,20);
    s=0;
    figure(i)
    while nrmv > 1e-3                                                            % L-scheme iterations
        A.op  = @(u) L*(u-u0)-(diff(u,2)-c*diff(u0)*(diff(u)+1)+exp(c*u0)*F/a);  % Stabilized linear operator - implicit L-scheme
        A.lbc = @(u) a*(diff(u)+1)-q*exp(c*u0);                                  % Neumann condition at x=x1
        A.rbc = @(u) u-b;                                                        % Dirichlect condition at x=x2
        u = A\0;                                                                 % Solve the linearized BVP
        v = u-u0; nrmv = norm(v);                                                % 2-norm of Newton update
        s=s+1;
        u0=u;
    end
    figure(100)
    plot(u,x,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5); hold on
    set(gca,'Ydir','reverse');set(gca,'XAxisLocation','top');
    xlabel('$\psi(z)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
    legend(leg,'Location','best')
    u_init=u;
    [degree]=Function_Chebop_Gardner_Fig2(i,u_init);
end
figure(101); legend(leg,'Box','off','Location','best')
toc

function [F]=step_h(x,h,q1)
F_neg = h*q1;
F_pos = 0;
neg = x <= h;
F = F_neg.*neg + F_pos.*~neg;
end

%% function for L-scheme/chebop switching

function [degree]=Function_Chebop_Gardner_Fig2(i,u_init)
cheboppref.setDefaults('plotting','off','display','iter','maxIter',100)
z1=0; z2=2; D=[z1 z2]; % domain
b=1;                   % constant head \psi=b on bottom
q=0.1;                 % constan flux on top
a=1; c=0.5;            % parameters for K(\psi)=a*exp(-c*\psi)
z = chebfun('z',D);
N = chebop(D);
ValueArray = {'o','+','x','*'}';
hh=[0 0.5 1 1.5 2];
figure(101); hold on

h=hh(i);
leg{i}=['h =',num2str(h)];
fprintf('h= %2.2f i= %d ***********************************************\n',h,i);
str=['h=',num2str(h)];
[F]=step_h(z,h,q);     % piecewise constant water uptake
N.op=@(x,u) diff(u,2)-c*diff(u).*(diff(u)+1)+exp(c*u)*F/a;
N.lbc=@(u) a*(diff(u)+1)-q*exp(c*u);
N.rbc=@(u) u-b;
N.init=u_init;
[u,info] = solvebvp(N,0);
plot(u,z,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',7);
set(gca,'Ydir','reverse');set(gca,'XAxisLocation','top');
xlabel('$\psi(z)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); xlim([1 inf]); ylim([0 2]);

op_residual = norm(N.op(u));                                  % Residual of the differential equation
lbc=N.lbc(u); lbc_residual =abs(lbc(z1));                     % Residuals of left_BC
rbc=N.rbc(u); rbc_residual =abs(rbc(z2));                     % Residuals of Right-BC
Residual = sqrt(op_residual^2+lbc_residual^2+rbc_residual^2); % Residual-BVP
degree=length(u);                                             % degree of Chebyshev polynomial
fprintf('op_residual = %2.2e \n',op_residual);
fprintf('lbc_residual = %2.2e  rbc_residual = %2.2e \n',lbc_residual,rbc_residual);
fprintf('Residual/norm(u) = %2.2e norm(u) = %2.2e\n',Residual/norm(u),norm(u));
fprintf('degree of Chebyshev polynomial= %d \n',degree);

xNewton=1:length(info.normDelta);
Newton_nq=info.normDelta(xNewton);
coeff=abs(chebcoeffs(u,degree));
figure(i); hold on
title(leg(i));
semilogy(xNewton,Newton_nq,'m*-','LineWidth',1,'MarkerSize',7);
xlabel('Iteration count $k$','Interpreter','latex');
ylabel('$\|\psi^{k}-\psi^{k-1}\|$','Interpreter','latex'); set(gca,'yscale','log');
figure(i+4); hold on
title(leg(i));
semilogy(coeff(1:1:end),'o-','LineWidth',1,'MarkerSize',7);
xlabel('Degree $j$ of the Chebyshev polinomial','Interpreter','latex');
ylabel('Magnitude $|\psi_{j}|$ of the coefficients','Interpreter','latex'); set(gca,'yscale','log');
cheboppref.setDefaults('factory'); % default
end

%% Elapsed time is 12.448614 seconds.
