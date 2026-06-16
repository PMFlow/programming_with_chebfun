%% Newton & explicit L-scheme - fixed-point problem f(x)=0;

close all; 
tic
iL=0; % Newton
% iL=1; % L-scheme
L=5;
x_ref = 0; a=10; 
f = @(x) x + x^2 + a*x^3;              % nonlinear problem f(x)=0;
x=0.5;                                 % Initial guess
nrmv = 1;                              % Initialize variables
corr=zeros(1,20); rez_v=zeros(1,20); rez_u=zeros(1,20); xx=zeros(1,20);
rez_uinit=zeros(1,20); err_u=zeros(1,20);
s=0;
while nrmv > 1e-6                      % Newton iterations
    d=1+2*x+3*a*x.^2;                  % d:=f'
    if iL==0
        v=-f(x)/d;        % <=== Solve v=-f(x)/f'(x); <=== f(x)+f'(x)*v=0;
    else
        v=-sign(d)*f(x)/L;             % <=== Solve f(x_k)=L*(x_{k+1}-x_k); (x_{k+1}-x_k)=v;
    end
    nrmv = norm(v);                    % 2-norm of Newton update
    x = x+v;
    if mod(s,10)==0; fprintf('s = %d  v= %2.2e  x= %2.2e  \n',s,v,x); end
    s=s+1; xx(s)=x;
    corr(s)=nrmv;
    err_u(s)=norm(x-x_ref);
    rez_u(s)=norm(f(x));
end
figure; 
dk=1; 
plot(xx(1:dk:s),'r*-')
xlabel('Iteration count $k$','Interpreter','latex');
ylabel('$x(k)$','Interpreter','latex');
figure;
semilogy(corr(1:dk:end),'m*-','LineWidth',1); hold on
semilogy(rez_u(1:dk:end),'gs-','LineWidth',1);
semilogy(err_u(1:dk:end),'b+-','LineWidth',1);
xlabel('Iteration count $k$','Interpreter','latex');
legend('Corrections = $\|v\|$','Residuals = $\|f(x)\|$','Errors = $\|x-x*\|$','Interpreter','latex');
legend('Box','off','Location','best')
toc
