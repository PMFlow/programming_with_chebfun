%% Lscheme_chebop: provides an initial guess by the implicit L-scheme and solves the BVP with 'chebop'
close all; clear all;
tic
D=[0 5];              % domain
A = chebop(D);        % linear operator
x = chebfun('x',D);   % independent variable
u0=chebfun(1,D);      % initial iterate 
nrmv = 1; rez=1;
s=1;
L=5;
phi = @(u) diff(u,2)+2*u.*sin(u);
S=20; corr=zeros(1,S); rez_u=zeros(1,S);
while nrmv > 1e-3 
    A.op  = @(u) L*(u-u0)-(diff(u,2)+2*u.*sin(u0));
    A.lbc = 'neumann';
    A.rbc = @(u) diff(u).*u0-2;
    u = A\0;
    v = u-u0; 
    nrmv = norm(v);   % norm of corrections  
    corr(s)=nrmv;
    rez=norm(phi(u)); % norm of residuals
    rez_u(s)=rez;
    if mod(s,1)==0
        fprintf('s = %d ||v||= %2.2e ||phi(u)|| %2.2e \n',s,nrmv,rez);
    end
    u0=u;
    s=s+1;
end
figure
semilogy(corr,'m*-','LineWidth',1); hold on
semilogy(rez_u,'b+-','LineWidth',1);
xlabel('Iteration count $k$','Interpreter','latex');
legend('Corrections = $||v||$','Residuals = $||\Phi(u)||$','Interpreter','latex');
legend('Box','off','Location','best')

% 
cheboppref.setDefaults('plotting','off','display','iter') % iteration details
N = chebop(D);                        % Create a chebop on [D]
N.op = @(u) diff(u,2)+2*u.*sin(u);    % Assign the DE part
N.lbc = 'neumann';                    % Homogeneous Neumann LBC
N.rbc = @(u) diff(u).*u-2;            % Set RBC, u'(5)*u(5) = 2
N.init = u0;
[u,info] = solvebvp(N,0);
degree=length(u); 
s=1:length(info.normDelta);
Newton=info.normDelta(s);
coeff=abs(chebcoeffs(u,degree));
figure;
semilogy(s,Newton,'m*-','LineWidth',1,'MarkerSize',7);
yticks([1e-15 1e-10 1e-5 1]);
xlabel('Iteration count $k$','Interpreter','latex');
ylabel('$\|u^{s}-u^{s-1}\|$','Interpreter','latex');
figure;
semilogy(coeff,'o-','LineWidth',1,'MarkerSize',7);
xlabel('Degree of the Chebyshev polinomial','Interpreter','latex');
ylabel('Amplitude of the coefficients','Interpreter','latex');
figure
plot(x,u,'ro-','LineWidth',1); hold on
xlabel('$x$','Interpreter','latex');
ylabel('$u(x)$','Interpreter','latex');
toc
