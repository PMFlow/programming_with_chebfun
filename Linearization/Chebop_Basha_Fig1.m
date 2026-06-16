%% Chebfun solutions with Basha model for different (n,q) - [Basha, 1999, Fig1]
clear all; close all; 
tic
cheboppref.setDefaults('plotting','off','display','iter','maxIter',200) % iteration details
ValueArray = {'o','+','x','*'}';
z1=0; z2=2; D=[z1 z2];         % solutions's domain
b=1;                           % constant head \psi=b on bottom boundary
z = chebfun('z',D);            % vertical coordinate as a 'chebfun'
N = chebop(D);
zz=0:0.02:2;
uu=zeros(4,length(zz));
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
    fprintf('n= %2.2f  q= %2.2f  i= %d ***********************************************\n',n,q,i);
    str=['n=',num2str(n), 'q=',num2str(q)];
    F= 0;
    N.op = @(x,u) (1+u^n)*diff(u,2)-n*(u^(n-1))*diff(u)*(diff(u)+1)+(1+u^n)^2*F;   % op. a=1 % changed wrong u^(n+1) into u^(n-1)!!!
    N.lbc=@(u) diff(u)-(1+u^n)*q+1;
    N.rbc=@(u) u-b;
    [u,info] = solvebvp(N,0);  % Solve using overloaded \
    plot(u,z,'LineWidth',1,'Marker',ValueArray{i},'MarkerSize',5);
    hold on
    leg{i}=([' n=',num2str(n) '   q=',num2str(q)]); 
    legend(leg,'Box','off','Location','best')
    op_residual = norm(N.op(u));              % Residual of the differential equation
    lbc=N.lbc(u); lbc_residual =abs(lbc(z1)); % Residuals of boundary conditions
    rbc=N.rbc(u); rbc_residual =abs(rbc(z2)); % Residuals of boundary conditions
    Residual = sqrt(op_residual^2+lbc_residual^2+rbc_residual^2); % Residual of the BVP
    degree=length(u);                         % degree of Chebyshev polynomial
    uu(i,:)=u(zz);
    fprintf('op_residual = %2.2e lbc_residual = %2.2e  rbc_residual = %2.2e \n', op_residual, lbc_residual, rbc_residual);
    fprintf('Residual = %2.2e Residual/norm(u) = %2.2e norm(u) = %2.2e\n',Residual,Residual/norm(u),norm(u));
    fprintf('degree of Chebyshev polynomial= %d \n',degree);
    Q=q*100;
end
hold off
ylabel('$z$','Interpreter','latex'),xlabel('$\psi(z)$','Interpreter','latex')
set(gca,'Ydir','reverse'); set(gca,'XAxisLocation','top')
xNewton=1:length(info.normDelta);
Newton_nq=info.normDelta(xNewton);
coeff=abs(chebcoeffs(u,degree));%
figure;
semilogy(xNewton,Newton_nq,'m*-','LineWidth',1,'MarkerSize',7);
yticks([1e-15 1e-10 1e-5 1]);
xlabel('Iteration count $k$','Interpreter','latex');
ylabel('$\|\psi^{k}-\psi^{k-1}\|$','Interpreter','latex');
figure;
semilogy(coeff,'o-','LineWidth',1,'MarkerSize',7);
xlabel('Degree $j$ of the Chebyshev polinomial','Interpreter','latex');
ylabel('Magnitude $|\psi_{j}|$ of the coefficients','Interpreter','latex');
toc
%% Elapsed time is 9.742923 seconds.
