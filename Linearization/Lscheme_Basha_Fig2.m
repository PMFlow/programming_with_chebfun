%% L-scheme solutions with Basha model for inceeasing uptake depth h

close all; clear all;  
% cheboppref.setDefaults('plotting','off','display','iter','maxIter',500) % iteration details
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
L=1e2; 
for i=1:4
    h=hh(i);
    fprintf('h= %2.2f i= %d ***********************************************\n',h,i);
    [F]=step_h(x,h,q);                        % piecewise constant water uptake 
    leg{i}=['h =',num2str(h)];
    phi = @(u) (1+u^n).*diff(u,2)-n*u.^(n-1).*diff(u).*(diff(u)+1)+(1+u.^n).^2.*F;
    u0=[chebfun(1,D)]; u = u0; u_init=u;      % Initial guess
    nrmv = 1; y = 0;                          % Initialize variables
    s=1; 
    figure(i); 
    plot3(x,chebfun(y,D),u,'LineWidth',0.75), hold on
    S=20; corr=zeros(1,S); rez_v=zeros(1,S); rez_u=zeros(1,S);
    while nrmv > 1e-11                                                                             % L-scheme iterations
        A.op = @(u) L*(u-u0)-((1+u0^n)*diff(u,2)-n*(u0^(n-1))*diff(u0)*(diff(u)+1)+(1+u0^n)^2*F);  % Stabilized linear operator
        A.lbc = @(u) diff(u)-(1+u0^n)*q+1;                                                         % Neumann condition at x=x1
        A.rbc = @(u) u-b;                                                                          % Dirichlect condition at x=x2
        u = A\0;                                                                                   % Solve the linearized BVP
        v = u-u0; nrmv = norm(v); y = y+norm(v);                                                   % Corrections norm
        corr(s)=nrmv;
        rez=norm(phi(u));
        rez_u(s)=rez;
        if mod(s,10)==0
            fprintf('s = %d    ||v||= %2.2e ||phi(u)|| %2.2e \n',s,nrmv,rez);
            plot3(x,chebfun(y,D),u)
            xlabel('$x$','Interpreter','latex');
            ylabel('$y=\sum{\|v_s\|}$','Interpreter','latex');
            zlabel('$u_s$','Interpreter','latex'); box on
        end
        u0=u;
        s=s+1;
    end
    figure(5+i);
    semilogy(corr,'m*-','LineWidth',1); hold on
    semilogy(rez_u,'b+-','LineWidth',1);
    title(leg(i));
    xlabel('Iteration count $s$','Interpreter','latex');
    legend('Corrections = $||v||$','Residuals = $||\Phi(u)||$','Interpreter','latex');
    legend('Box','off','Location','best')
    figure(112)
    plot(u,x,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5); hold on
end
hold off
legend(leg,'Box','off','Location','best')
ylabel('$z$','Interpreter','latex'),xlabel('$\psi(z)$','Interpreter','latex')
set(gca,'Ydir','reverse'); set(gca,'XAxisLocation','top'); ylim([0 2]);
toc
function [F]=step_h(x,h,q1)
F_neg = h*q1;
F_pos = 0;
neg = x <= h;
F = F_neg.*neg + F_pos.*~neg;
end
% Elapsed time is 79.808779 seconds.
