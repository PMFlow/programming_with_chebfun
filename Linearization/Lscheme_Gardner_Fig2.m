%% L-scheme for inhomogeneous BVP with Grdner model (Eqs. 16 and 17)


tic
clear all; close all; % clc;
% cheboppref.setDefaults('plotting','off','display','iter','maxIter',500)
cheboppref.setDefaults('factory'); % default
x1=0; x2=2; D=[x1 x2]; xx=0:0.02:2; L=0.5;
b=1;      % constant head \psi=b on bottom
q=0.1;    % constant flux on top
a=1; c=0.5; % parameters for K(\psi)=a*exp(-c*\psi)
A = chebop(D);
x = chebfun('x',D);
ValueArray = {'o','+','x','*','s'}';
hh=[0 0.5 1 1.5];

for i=1:4
    h=hh(i);
    fprintf('h= %2.2f i= %d ***********************************************\n',h,i);
    leg{i}=['h =',num2str(h)];
    [F]=step_h(x,h,q);                                        % piecewise constant water uptake
    phi = @(u) diff(u,2)-c*diff(u).*(diff(u)+1)+exp(c*u)*F/a; % ODE part of the nonlinear BVP
    % u=chebfun(1.,D); u_init=u; u0=u;                      % Initial guess
    u = x; u_init=u;  u0=u;                                   % Initial guess
    nrmv = 1; y = 0;                                          % Initialize variables
    corr=zeros(1,20); rez_v=zeros(1,20); rez_u=zeros(1,20); rez_uinit=zeros(1,20);
    s=0;
    figure(i)
    plot3(x,chebfun(y,D),u,'LineWidth',0.75), hold on
    while nrmv > 1e-11                                                           % Newton iterations
        A.op  = @(u) L*(u-u0)-(diff(u,2)-c*diff(u0)*(diff(u)+1)+exp(c*u0)*F/a);  % Stabilized linear operator - implicit L-scheme
        A.lbc = @(u) a*(diff(u)+1)-q*exp(c*u0);                                  % Neumann condition at x=x1
        A.rbc = @(u) u-b;                                                        % Dirichlect condition at x=x2
        u = A\0;                                                                 % Solve the linearized BVP
        v = u-u0; nrmv = norm(v); y = y+norm(v);                                 % 2-norm of Newton update
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
    figure(5+i);
    semilogy(corr,'m*-','LineWidth',1); hold on
    semilogy(rez_u,'b+-','LineWidth',1);
    title(leg(i));
    xlabel('Iteration count $s$','Interpreter','latex');
    legend('Corrections = $||v||$','Residuals = $||\Phi(u)||$','Interpreter','latex');
    legend('Box','off','Location','southwest')
    figure(112)
    plot(u,x,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5); hold on
end
legend(leg,'Location','best')
set(gca,'yscale','linear');set(gca,'Ydir','reverse');set(gca,'XAxisLocation','top');
xlabel('$\psi(z)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
toc
function [F]=step_h(x,h,q1)
F_neg = h*q1;
F_pos = 0;
neg = x <= h;
F = F_neg.*neg + F_pos.*~neg;
end
% Elapsed time is 21.617839 seconds.
