%% Newton-Frechet solutions with Basha model for different (n,q) - [Basha, 1999, Fig1]

close all; clear all;
cheboppref.setDefaults('plotting','off','display','iter','maxIter',500) % iteration details 
% cheboppref.setDefaults('factory')
tic
x1=0; x2=2; D=[x1 x2]; xx=0:0.02:2; L=0.5;
A = chebop(D);           % Operator on the interval [0,5]
x = chebfun('x',D);      % The chebfun "x" on the interval [0,5]
ValueArray = {'o','+','x','*'}';
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
    fprintf('n= %2.2f  q= %2.2f  i= %d ***************************************\n',n,q,i);,'LineWidth',0.75
    leg{i}=['n=',num2str(n), '   q=',num2str(q)];
    if n==3
        phi = @(u) (1+u.^3).*diff(u,2)-3*u.^2.*diff(u).*(diff(u)+1); % ODE part of the nonlinear BVP - n=3
    else
        phi = @(u) (1+u.^7).*diff(u,2)-7*u.^6.*diff(u).*(diff(u)+1); % ODE part of the nonlinear BVP - n=7
    end
    u=chebfun(0.,D); u_init=u; % !!! OK !!!
    nrmv = 1; y = 0;                            % Initialize variables
    s=0;
    figure(i);
    plot3(x,chebfun(y,D),u,'LineWidth',0.75), hold on
    S=20; corr=zeros(1,S); rez_v=zeros(1,S); rez_u=zeros(1,S);
    while nrmv > 1e-11                          % Newton iterations
        fprintf('Iteration count= %d \n',s+1);
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
        nrmv = norm(v); y = y+norm(v);                                                                                   % 2-norm of Newton update
        u = u+v;
        s=s+1;
        corr(s)=nrmv;
        rez_v(s)=norm(A.op(v));
        rez_u(s)=norm(phi(u));
        plot3(x,chebfun(y,D),u)
        xlabel('$x$','Interpreter','latex');
        ylabel('$y=\sum{\|v_s\|}$','Interpreter','latex');
        zlabel('$u_s$','Interpreter','latex'); box on
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
set(gca,'Ydir','reverse'); set(gca,'XAxisLocation','top')
toc
% Elapsed time is 11.705256 seconds.
