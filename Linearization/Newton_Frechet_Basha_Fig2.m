%% Newton-Frechet solutions with Basha model for inceeasing uptake depth h

close all; clear all;  % clc;
% cheboppref.setDefaults('plotting','off','display','iter','maxIter',500) % 200) % iteration details - adaugat !!!
cheboppref.setDefaults('factory')
tic
n=3;                             % exponent in the constitutive law K(\psi)=1/(1+\psi^n)
q=0.1;                           % constant flux on top boundary
b=1;                             % constant head \psi=b on bottom boundary
x1=0; x2=2; D=[x1 x2];
A = chebop(D);                   % Operator on the interval [0,5]
x = chebfun('x',D);              % The chebfun "x" on the interval [0,5]
ValueArray = {'o','+','x','*','s'}';
hh=[0 0.5 1 1.5]; 
for i=1:4
    h=hh(i);
    fprintf('h= %2.2f i= %d ***********************************************\n',h,i);
    str=['h=',num2str(h)];
    [F]=step_h(x,h,q);                                                              % piecewise constant water uptake 
    leg{i}=['h =',num2str(h)];
    if n==2
        phi = @(u) (1+u.^2).*diff(u,2)-2*u.*diff(u).*(diff(u)+1)+(1+u.^2).^2.*F;    % ODE part of the nonlinear BVP - n=2
    elseif n==3
        phi = @(u) (1+u.^3).*diff(u,2)-3*u.^2.*diff(u).*(diff(u)+1)+(1+u.^3).^2.*F; % ODE part of the nonlinear BVP - n=3
    elseif n==5
        phi = @(u) (1+u.^5).*diff(u,2)-5*u.^4.*diff(u).*(diff(u)+1)+(1+u.^5).^2.*F; % ODE part of the nonlinear BVP - n=5
    elseif n==7
        phi = @(u) (1+u.^7).*diff(u,2)-7*u.^6.*diff(u).*(diff(u)+1)+(1+u.^7).^2.*F; % ODE part of the nonlinear BVP - n=7
    end
    u=chebfun(0.,D);  
    nrmv = 1; y = 0;  rez=1;  
    s=0; 
    figure(i);
    plot3(x,chebfun(y,D),u,'LineWidth',0.75), box on; hold on    
    title(leg(i));
    S=20; corr=zeros(1,S); rez_v=zeros(1,S); rez_u=zeros(1,S);
    while nrmv > 1e-11 
        Du = diff(u); 
        if n==2
            A.op = @(v) (1+u^2)*diff(v,2)-(4*u*diff(u)+2*u)*diff(v)+(2*u*diff(u,2)+2*diff(u)*(diff(u)+1)+(2*u+4*u^3)*F)*v;
            A.lbc = @(v) diff(v)-2*u*q*v+Du(x1)-(1+u(x1)^2)*q+1;
        elseif n==3
            A.op = @(v) (1+u^3)*diff(v,2)-(6*u^2*diff(u)+3*u^2)*diff(v)+(3*u^2*diff(u,2)-6*u*diff(u)*(diff(u)+1)+3*u^2*F)*v;     % Frechet derivative at u - n=3
            A.lbc = @(v) diff(v)-3*u(x1)^2*v*q+Du(x1)-(1+u(x1)^3)*q+1;                                                           % Neumann condition at x=x1 - n=3
        elseif n==5
            A.op = @(v) (1+u^5)*diff(v,2)-(10*u^4*diff(u)+5*u^4)*diff(v)+(5*u^4*diff(u,2)-20*u^3*diff(u)*(diff(u)+1)+5*u^4*F)*v; % Frechet derivative - n=5
            A.lbc = @(v) diff(v)-5*u(x1)^4*q*v+Du(x1)-(1+u(x1)^5)*q+1;                                                           % Neumann condition at x=x1 - n=5
        elseif n==7
            A.op = @(v) (1+u^7)*diff(v,2)-(14*u^6*diff(u)+7*u^6)*diff(v)+(7*u^6*diff(u,2)-42*u^5*diff(u)*(diff(u)+1)+7*u^6*F)*v; % Frechet derivative at u - n=7
            A.lbc = @(v) diff(v)-7*u(x1)^6*v*q+Du(x1)-(1+u(x1)^7)*q+1;                                                           % Neumann condition at x=x1 - n=7
        end
        A.rbc = @(v) v+u(x2)-1;                                                                                                  % Dirichlect condition at x=x2        
        A.rbc = @(v) v+u(x2)-1;         % Dirichlect condition at x=x2
        v = A\(-phi(u));                % Solve the linearized BVP
        nrmv = norm(v); y = y+norm(v);  % 2-norm of Newton update
        u = u+v;
        plot3(x,chebfun(y,D),u); hold on
        s=s+1; 
        corr(s)=nrmv;
        rez_v(s)=norm(A.op(v)); 
        rez=norm(phi(u));
        rez_u(s)=rez; 
        if mod(s,10)==0
            fprintf('s = %d    ||v||= %2.2e ||phi(u)|| %2.2e \n',s,nrmv,rez);
            plot3(x,chebfun(y,D),u)
            xlabel('$x$','Interpreter','latex');
            ylabel('$y=\sum{\|v_s\|}$','Interpreter','latex');
            zlabel('$u_s$','Interpreter','latex'); box on            
        end
    end
    figure(4+i);
    semilogy(corr,'m*-','LineWidth',1); hold on
    semilogy(rez_u,'b+-','LineWidth',1);
    semilogy(rez_v,'gx-','LineWidth',1);
    title(leg(i));
    xlabel('Iteration count $s$','Interpreter','latex');
    legend('Corrections = $||v||$','Residuals u = $||\Phi(u)||$','Residuals v = $||A.op(v)||$','Interpreter','latex');
    legend('Box','off','Location','best')
    figure(111)
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
% Elapsed time is 49.316568 seconds.
