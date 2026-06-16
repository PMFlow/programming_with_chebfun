%% L-scheme solutions with Basha model for different (n,q) - [Basha, 1999, Fig1]

close all; clear all;  % clc;
% cheboppref.setDefaults('plotting','off','display','iter','maxIter',500) % iteration details 
cheboppref.setDefaults('factory')
tic
b=1;                        % constant head \psi=b on bottom
x1=0; x2=2; D=[x1 x2]; xx=0:0.02:2;
A = chebop(D);              % Operator on the interval [0,5]
x = chebfun('x',D);         % The chebfun "x" on the problem interval
ValueArray = {'o','+','x','*'}';
L=100;
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
    fprintf('n= %2.2f  q= %2.2f  i= %d ***************************************\n',n,q,i);
    leg{i}=['n=',num2str(n), '   q=',num2str(q)];
    phi = @(u) (1+u^n)*diff(u,2)-n*(u^(n-1))*diff(u)*(diff(u)+1);                     % ODE part of the nonlinear BVP
    u0=chebfun(1,D); u = u0;                                                        % Initial guess
    nrmv = 1; y = 0;  rez=1;                                                          % Initialize variables
    s=1; 
    figure(i)
    plot3(x,chebfun(y,D),u,'LineWidth',0.75), hold on
    S=20; corr=zeros(1,S); rez_v=zeros(1,S); rez_u=zeros(1,S);
    while nrmv > 1e-11
        A.op  = @(u) L*(u-u0)-((1+u0^n)*diff(u,2)-n*(u0^(n-1))*diff(u0)*(diff(u)+1));  % Stabilized linear operator
        A.lbc = @(u) diff(u)-(1+u0^n)*q+1;                                             % Neumann condition at x=x1
        A.rbc = @(u) u-b;                                                              % Dirichlect condition at x=x2
        u = A\0;                                                                       % Solve the linearized BVP
        v = u-u0; nrmv = norm(v); y = y+norm(v);                                       % 2-norm of Newton update
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
    figure(4+i);
    semilogy(corr,'m*-','LineWidth',1); hold on
    semilogy(rez_u,'b+-','LineWidth',1);
    title(leg(i));
    xlabel('Iteration count $s$','Interpreter','latex');
    legend('Corrections = $||v||$','Residuals = $||\Phi(u)||$','Interpreter','latex');
    legend('Box','off','Location','best')
    figure(111)
    plot(u,x,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5); hold on
end
hold off
legend(leg,'Box','off','Location','best')
ylabel('$z$','Interpreter','latex'),xlabel('$\psi(z)$','Interpreter','latex')
set(gca,'Ydir','reverse'); set(gca,'XAxisLocation','top')
toc
% Elapsed time is 44.715517 seconds.
