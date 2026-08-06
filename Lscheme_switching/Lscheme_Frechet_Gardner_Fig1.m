%% L-scheme/Frechet switching for homogeneous BVP with Grdner model (Eqs. 16 and 17)

tic
clear all; close all;
cheboppref.setDefaults('factory'); % default 
x1=0; x2=2; D=[x1 x2]; xx=0:0.02:2;
b=1;      % constant head \psi=b on bottom
q=0.1;    % constant flux on top
a=1; c=2; % parameters for K(\psi)=a*exp(-c*\psi)
A = chebop(D);
x = chebfun('x',D);
ValueArray = {'o','+','x','*'}';
L=2;
for i=1:4
    c=0.5*i;
    leg{i}=['c=',num2str(c)];
    phi = @(u) diff(u,2)-c*diff(u).*(diff(u)+1); % ODE part of the nonlinear BVP
    u = x; u0=u;                                 % Initial guess
    nrmv = 1;                        
    s=0;
    while nrmv > 1e-3                                              % L-scheme iterations
        A.op  = @(u) L*(u-u0)-(diff(u,2)-c*diff(u0)*(diff(u)+1));  % Stabilized linear operator - implicit L-scheme
        A.lbc = @(u) a*(diff(u)+1)-q*exp(c*u0);                    % Neumann condition at x=x1
        A.rbc = @(u) u-b;                                          % Dirichlect condition at x=x2
        u = A\0;                                                   % Solve the linearized BVP
        v = u-u0; nrmv = norm(v);                                  % 2-norm of Newton update
        s=s+1;
        u0=u;
    end
    figure(100)
    plot(u,x,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',7); hold on
    set(gca,'Ydir','reverse');set(gca,'XAxisLocation','top'); 
    xlabel('$\psi(z)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); ylim([0 2]);
    legend(leg,'box','off','Location','best')
    u_init=u;
    Function_Frechet_Gardner_Fig1(i,u_init);
end
figure(101); legend(leg,'Box','off','Location','best')
toc

%% function for L-scheme/Frechet switching

function Function_Frechet_Gardner_Fig1(i,u_init)

cheboppref.setDefaults('factory'); % default
x1=0; x2=2; D=[x1 x2]; xx=0:0.02:2; L=100;
b=1;      % constant head \psi=b on bottom
q=0.1;    % constant flux on top
a=1; c=2; % parameters for K(\psi)=a*exp(-c*\psi)
A = chebop(D);
x = chebfun('x',D);
ValueArray = {'o','+','x','*'}';
    c=0.5*i;
    fprintf('c= %2.2f  i= %d ***************************************\n',c,i);
    leg{i}=['c=',num2str(c)];
    phi = @(u) diff(u,2)-c*diff(u).*(diff(u)+1);      % ODE part of the nonlinear BVP
    u = u_init;                                       % Initial guess
    nrmv = 1;                                         
    corr=zeros(1,20); rez_v=zeros(1,20); rez_u=zeros(1,20); 
    s=0;
    while nrmv > 1e-11                                                           % Newton iterations
        A.op = @(v) diff(v,2)-c*(2*diff(u)+1)*diff(v);                           % Frechet differential at u
        Du = diff(u);                                                            % Needed to compute the BCs
        A.lbc = @(v) a*diff(v)-q*exp(c*u(x1))*(c*v)+a*(Du(x1)+1)-q*exp(c*u(x1)); % Neumann condition at x=x1
        A.rbc = @(v) v+u(x2)-1;                                                  % Dirichlect condition at x=x2
        v = A\(-phi(u));                                                         % Solve the linearized BVP
        nrmv = norm(v);                                                          % 2-norm of Newton update
        u = u+v;
        title(leg(i));
        xlabel('$x$','Interpreter','latex');
        ylabel('$y=\sum{\|v_s\|}$','Interpreter','latex');
        zlabel('$u_s$','Interpreter','latex'); box on
        s=s+1;
        corr(s)=nrmv;
        rez_v(s)=norm(A.op(v));
        rez_u(s)=norm(phi(u));
    end
    figure(i);
    semilogy(corr,'m*-','LineWidth',1); hold on
    semilogy(rez_u,'b+-','LineWidth',1);
    semilogy(rez_v,'gs-','LineWidth',1);
    title(leg(i));
    xlabel('Iteration count $k$','Interpreter','latex');
    legend('Corrections = $||v||$','Residuals u = $||\Phi(u)||$','Residuals v = $||A.op(v)||$','Interpreter','latex');
    legend('Box','off','Location','best')
    figure(101)
    if i==1
        plot(u,x,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5); hold on
    else
        plot(u(xx(1:4:end)),xx(1:4:end),'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5); hold on
    end
set(gca,'yscale','linear');set(gca,'Ydir','reverse');set(gca,'XAxisLocation','top');
xlabel('$\psi(z)$','Interpreter','latex'); ylabel('$z$','Interpreter','latex'); ylim([0 2]);
fprintf('op_residual = %2.2e \n',rez_u(s));
end

%% Elapsed time is 11.638095 seconds.
