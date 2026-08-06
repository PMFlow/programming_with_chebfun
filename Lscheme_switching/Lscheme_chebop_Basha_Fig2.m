%% L-scheme/chebop switching solutions with Basha model for inceeasing uptake depth h  (Eqs. 16 and 18)

close all; clear all;
cheboppref.setDefaults('factory') 
tic
n=3;                                    % exponent in the constitutive law K(\psi)=1/(1+\psi^n)
q=0.1;                                  % constant flux on top boundary
b=1;                                    % constant head \psi=b on bottom boundary
x1=0; x2=2; D=[x1 x2];
A = chebop(D);                          % Operator on the interval [0,5]
x = chebfun('x',D);                     % The chebfun "x" on the problem interval
ValueArray = {'o','+','x','*','s'}';
hh=[0 0.5 1 1.5 2];
L=100;
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
    xlabel('$\psi(z)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
    legend(leg,'Box','off','Location','best')
    u_init=u;
    [degree]=Function_Chebop_Basha_Fig2(i,u_init);
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

function [degree]=Function_Chebop_Basha_Fig2(i,u_init)
cheboppref.setDefaults('plotting','off','display','iter','maxIter',500) % iteration details
z1=0; z2=2; D=[z1 z2]; % solutions's domain
n=3;                   % exponent in the constitutive law K(\psi)=1/(1+\psi^n)
q=0.1;                 % constant flux on top boundary
b=1;                   % constant head \psi=b on bottom boundary
z = chebfun('z',D);    % vertical coordinate as a 'chebfun'
N = chebop(D);
hh=[0 0.5 1 1.5 2];
ValueArray = {'o','+','x','*'}';
figure(101); hold on

h=hh(i);
leg{i}=(['h = ', num2str(h)]);
fprintf('h= %2.2f i= %d ***********************************************\n',h,i);
[F]=step_h(z,h,q);                                                % piecewise constant water uptake
N.op = @(z,u) (1+u^n)*diff(u,2)-n*(u^(n-1))*diff(u)*(diff(u)+1)+(1+u^n)^2*F;
N.lbc=@(u) diff(u)-(1+u^n)*q+1;
N.rbc=@(u) u-b;
N.init=u_init;
[u,info] = solvebvp(N,0);
plot(u,z,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5);
set(gca,'Ydir','reverse');set(gca,'XAxisLocation','top');
xlabel('$\psi(z)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); xlim([1 inf]); ylim([0 2]);

op_residual = norm(N.op(u));                                      % Residual of the differential equation
lbc=N.lbc(u); lbc_residual =abs(lbc(z1));                         % Residuals of boundary conditions
rbc=N.rbc(u); rbc_residual =abs(rbc(z2));                         % Residuals of boundary conditions
Residual = sqrt(op_residual^2+lbc_residual^2+rbc_residual^2);     % Residual of the BVP
degree=length(u);                                                 % degree of Chebyshev polynomial
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
xlabel('Iteration count $s$','Interpreter','latex');
ylabel('$\|\psi^{k}-\psi^{k-1}\|$','Interpreter','latex'); set(gca,'yscale','log');
figure(i+4); hold on
title(leg(i));
semilogy(coeff,'o-','LineWidth',1,'MarkerSize',7);
xlabel('Degree $j$ of the Chebyshev polinomial','Interpreter','latex');
ylabel('Magnitude $|\psi_{j}|$ of the coefficients','Interpreter','latex'); set(gca,'yscale','log');
cheboppref.setDefaults('factory'); % default
end
%% Elapsed time is 31.375830 seconds.