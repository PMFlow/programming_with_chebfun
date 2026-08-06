%% L-scheme/chebop switching solutions with Basha model for different (n,q) - [Basha, 1999, Fig1]  (Eqs. 16 and 18)

close all; clear all;  
cheboppref.setDefaults('factory') 
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
    nrmv = 1; y = 0;  rez=1;                                                           
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
    [degree]=Function_Chebop_Basha_Fig1(i,u_init);
end
figure(101); legend(leg,'Box','off','Location','best')
toc

%% function for L-scheme/chebop switching

function [degree]=Function_Chebop_Basha_Fig1(i,u_init)
cheboppref.setDefaults('plotting','off','display','iter','maxIter',200) % iteration details
z1=0; z2=2; D=[z1 z2];         % solutions's domain
b=1;                           % constant head \psi=b on bottom boundary
z = chebfun('z',D);            % vertical coordinate as a 'chebfun'
N = chebop(D);
if i== 1
    n=3; q=0.01;
elseif i==2
    n=3; q=0.1;
elseif i==3
    n=7; q=0.01;
else
    n=7; q=0.1;
end
ValueArray = {'o','+','x','*'}';
figure(101); hold on

leg{i}=([' n=',num2str(n) '   q=',num2str(q)]);
fprintf('n= %2.2f  q= %2.2f  i= %d ***********************************************\n',n,q,i);
F= 0;
N.op = @(x,u) (1+u^n)*diff(u,2)-n*(u^(n-1))*diff(u)*(diff(u)+1)+(1+u^n)^2*F; 
N.lbc=@(u) diff(u)-(1+u^n)*q+1;
N.rbc=@(u) u-b;
N.init=u_init;
[u,info] = solvebvp(N,0);   
plot(u,z,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',7);
set(gca,'Ydir','reverse');set(gca,'XAxisLocation','top');
xlabel('$\psi(z)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); xlim([1 inf]); ylim([0 2]);

op_residual = norm(N.op(u));              % Residual of the differential equation
lbc=N.lbc(u); lbc_residual =abs(lbc(z1)); % Residuals of boundary conditions
rbc=N.rbc(u); rbc_residual =abs(rbc(z2)); % Residuals of boundary conditions
Residual = sqrt(op_residual^2+lbc_residual^2+rbc_residual^2); % Residual of the BVP
degree=length(u);                         % degree of Chebyshev polynomial
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
%% Elapsed time is 23.245111 seconds.