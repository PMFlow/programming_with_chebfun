%% 'chebop' implementation for inhomogeneous BVP with Grdner model (Eqs. 16 and 17)

tic
clear all; close all; 
cheboppref.setDefaults('plotting','off','display','iter','maxIter',500)
z1=0; z2=2; D=[z1 z2]; % domain
b=1;                   % constant head \psi=b on bottom
q=0.1;                 % constan flux on top
a=1; c=0.5;            % parameters for K(\psi)=a*exp(-c*\psi)
z = chebfun('z',D);
N = chebop(D);
ValueArray = {'o','+','x','*'}';
hh=[0 0.5 1 1.5 2];
figure;
for i=1:4
    h=hh(i);
    fprintf('h= %2.2f i= %d ***********************************************\n',h,i);
    str=['h=',num2str(h)];
    [F]=step_h(z,h,q);                     % piecewise constant water uptake
    N.op=@(x,u) diff(u,2)-c*diff(u).*(diff(u)+1)+exp(c*u)*F/a;
    N.lbc=@(u) a*(diff(u)+1)-q*exp(c*u);
    N.rbc=@(u) u-b;
    N.init=chebfun(0,D);                   % convergence
    % N.init=chebfun(1,D);                 % convergence
    % N.init=z;                            % no convergence
    [u,info] = solvebvp(N,0);
    plot(u,z,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5);
    hold on
    leg{i}=(['h = ', num2str(h)]);
    legend(leg,'Box','off','Location','best')
    op_residual = norm(N.op(u));                                  % Residual of the differential equation
    lbc=N.lbc(u); lbc_residual =abs(lbc(z1));                     % Residuals of left_BC
    rbc=N.rbc(u); rbc_residual =abs(rbc(z2));                     % Residuals of Right-BC
    Residual = sqrt(op_residual^2+lbc_residual^2+rbc_residual^2); % Residual-BVP
    degree=length(u);                                             % degree of Chebyshev polynomial
    fprintf('op_residual = %2.2e \n',op_residual);
    fprintf('lbc_residual = %2.2e  rbc_residual = %2.2e \n',lbc_residual,rbc_residual);
    fprintf('Residual/norm(u) = %2.2e norm(u) = %2.2e\n',Residual/norm(u),norm(u));
    fprintf('degree of Chebyshev polynomial= %d \n',degree);
end
hold off
ylabel('$z$','Interpreter','latex'),xlabel('$\psi(z)$','Interpreter','latex')
set(gca,'Ydir','reverse'); set(gca,'XAxisLocation','top'); 
xNewton=1:length(info.normDelta);
Newton_nq=info.normDelta(xNewton);
coeff=abs(chebcoeffs(u,degree));
figure;
semilogy(xNewton,Newton_nq,'m*-','LineWidth',1,'MarkerSize',7);
xlabel('Iteration count $k$','Interpreter','latex');
ylabel('$\|\psi^{k}-\psi^{k-1}\|$','Interpreter','latex');
figure;
semilogy(coeff(1:1:end),'o-','LineWidth',1,'MarkerSize',7);
xlabel('Degree $j$ of the Chebyshev polinomial','Interpreter','latex');
ylabel('Magnitude $|\psi_{j}|$ of the coefficients','Interpreter','latex');
toc
function [F]=step_h(x,h,q1)
F_neg = h*q1;
F_pos = 0;
neg = x <= h;
F = F_neg.*neg + F_pos.*~neg;
end
% Elapsed time is 6.709591 seconds for N.init=0
% no convergence for default N.init and N.init=z