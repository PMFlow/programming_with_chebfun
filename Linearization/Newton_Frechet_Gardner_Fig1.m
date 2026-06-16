%% Newton-Frechet method for homogeneous BVP with Grdner model (Eqs. 16 and 17)

tic
clear all; close all;
% cheboppref.setDefaults('plotting','off','display','iter','maxIter',500)
cheboppref.setDefaults('factory'); % default
x1=0; x2=2; D=[x1 x2]; xx=0:0.02:2; L=100;
b=1;      % constant head \psi=b on bottom
q=0.1;    % constant flux on top
a=1; c=2; % parameters for K(\psi)=a*exp(-c*\psi)
A = chebop(D);
x = chebfun('x',D);
ValueArray = {'o','+','x','*'}';
for i=1:4
    c=0.5*i;
    fprintf('c= %2.2f  i= %d ***************************************\n',c,i);
    leg{i}=['c=',num2str(c)];
    phi = @(u) diff(u,2)-c*diff(u).*(diff(u)+1); % ODE part of the nonlinear BVP
    u = x; u0=u;                                 % Initial guess
    nrmv = 1; y = 0;                             % Initialize variables
    corr=zeros(1,20); rez_v=zeros(1,20); rez_u=zeros(1,20); rez_uinit=zeros(1,20);
    s=0;
    figure(i)
    plot3(x,chebfun(y,D),u,'LineWidth',0.75), hold on
    while nrmv > 1e-11                                                           % Newton iterations
        A.op = @(v) diff(v,2)-c*(2*diff(u)+1)*diff(v);                           % Frechet differential at u
        Du = diff(u);                                                            % Needed to compute the BCs
        A.lbc = @(v) a*diff(v)-q*exp(c*u(x1))*(c*v)+a*(Du(x1)+1)-q*exp(c*u(x1)); % Neumann condition at x=x1
        A.rbc = @(v) v+u(x2)-1;                                                  % Dirichlect condition at x=x2
        v = A\(-phi(u));                                                         % Solve the linearized BVP
        nrmv = norm(v); y = y+norm(v);                                           % 2-norm of Newton update
        u = u+v;
        plot3(x,chebfun(y,D),u)
        title(leg(i));
        xlabel('$x$','Interpreter','latex');
        ylabel('$y=\sum{\|v_s\|}$','Interpreter','latex');
        zlabel('$u_s$','Interpreter','latex'); box on
        s=s+1;
        corr(s)=nrmv;
        rez_v(s)=norm(A.op(v));
        rez_u(s)=norm(phi(u));
        u0=u;
    end
    figure(4+i);
    semilogy(corr,'m*-','LineWidth',1); hold on
    semilogy(rez_u,'b+-','LineWidth',1);
    semilogy(rez_v,'gs-','LineWidth',1);
    title(leg(i));
    xlabel('Iteration count $k$','Interpreter','latex');
    legend('Corrections = $||v||$','Residuals u = $||\Phi(u)||$','Residuals v = $||A.op(v)||$','Interpreter','latex');
    legend('Box','off','Location','best')
    figure(112)
    if i==1
        plot(u,x,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5); hold on
    else
        plot(u(xx(1:4:end)),xx(1:4:end),'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5); hold on
    end
end
legend(leg,'Box','off','Location','best')
set(gca,'yscale','linear');set(gca,'Ydir','reverse');set(gca,'XAxisLocation','top');
xlabel('$\psi(z)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); ylim([0 2]);
toc
% Elapsed time is 17.033763 seconds.

