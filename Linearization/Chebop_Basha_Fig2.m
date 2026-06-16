%% Chebfun solutions with Basha model for inceeasing uptake depth h

clear all; close all; 
tic
cheboppref.setDefaults('plotting','off','display','iter','maxIter',500) % iteration details
% cheboppref.setDefaults('factory'); % default
ValueArray = {'o','+','x','*'}';
z1=0; z2=2; D=[z1 z2]; % solutions's domain
n=3;                   % exponent in the constitutive law K(\psi)=1/(1+\psi^n)
q=0.1;                 % constant flux on top boundary
b=1;                   % constant head \psi=b on bottom boundary
z = chebfun('z',D);    % vertical coordinate as a 'chebfun'
N = chebop(D);
zz=0:0.02:2;
uu=zeros(5,length(zz));
hh=[0 0.5 1 1.5 2];
figure;
for i=1:4
    h=hh(i);
    fprintf('h= %2.2f i= %d ***********************************************\n',h,i);
    [F]=step_h(z,h,q);                                                % piecewise constant water uptake 
    N.op = @(z,u) (1+u^n)*diff(u,2)-n*(u^(n-1))*diff(u)*(diff(u)+1)+(1+u^n)^2*F;
    N.lbc=@(u) diff(u)-(1+u^n)*q+1;
    N.rbc=@(u) u-b;
    N.init=chebfun(0.,D);  
    [u,info] = solvebvp(N,0);  
    u = simplify(u);
    plot(u,z,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5);
    hold on
    leg{i}=(['h = ', num2str(h)]); 
    legend(leg,'Location','best')
    op_residual = norm(N.op(u));                                      % Residual of the differential equation
    lbc=N.lbc(u); lbc_residual =abs(lbc(z1));                         % Residuals of boundary conditions
    rbc=N.rbc(u); rbc_residual =abs(rbc(z2));                         % Residuals of boundary conditions
    Residual = sqrt(op_residual^2+lbc_residual^2+rbc_residual^2);     % Residual of the BVP
    degree=length(u);                                                 % degree of Chebyshev polynomial
    fprintf('op_residual = %2.2e lbc_residual = %2.2e  rbc_residual = %2.2e \n', op_residual, lbc_residual, rbc_residual);
    fprintf('Residual = %2.2e Residual/norm(u) = %2.2e norm(u) = %2.2e\n',Residual,Residual/norm(u),norm(u));
    fprintf('degree of Chebyshev polynomial= %d \n',degree);
    uu(i,:)=u(zz);
end
hold off
ylabel('$z$','Interpreter','latex'),xlabel('$\psi(z)$','Interpreter','latex')
set(gca,'Ydir','reverse'); set(gca,'XAxisLocation','top')
xNewton=1:length(info.normDelta);
Newton_h=info.normDelta(xNewton);
coeff=abs(chebcoeffs(u,degree));
figure;
semilogy(xNewton,Newton_h,'m*-','LineWidth',1,'MarkerSize',7);
yticks([1e-15 1e-10 1e-5 1]);
xlabel('Iteration count $s$','Interpreter','latex');
ylabel('$\|\psi^{s}-\psi^{s-1}\|$','Interpreter','latex');
figure;
semilogy(coeff,'o-','LineWidth',1,'MarkerSize',7);
xlabel('Degree $j$ of the Chebyshev polinomial','Interpreter','latex');
ylabel('Magnitude $|\psi_{j}|$ of the coefficients','Interpreter','latex');
toc
function [F]=step_h(x,h,q1)
F_neg = h*q1;
F_pos = 0;
neg = x <= h;
F = F_neg.*neg + F_pos.*~neg;
end
% Elapsed time is 445.758802 seconds.

% h=0 and h=1: "Warning: Newton iteration failed. 
% Please try supplying a better initial guess via the .init field
% of the chebop. "

% h=1.5: "Warning: Newton iteration failed. Maximum number of iterations exceeded".
