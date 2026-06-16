%%   Estimated order of accuracy for NON-steady Richards-1D

close all; clear all;
tic
T=1;
y1=0; y2=1;
D=[0 T y1 y2];
x = chebfun2(@(x,y) x,D);
y = chebfun2(@(x,y) y,D);
E=chebfun2(@(x,y) -x.*y.*(1-y)-0.5,[0 1 0 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi0 = @(u) diffx(u)./((1-u).^2)-u.^2.*diffy(u,2)-2*u.*diffy(u).*(diffy(u)+1);
phis=phi0(E);      % source term as chebfun
Itest = 1;
S= 10000; iS=1:S;
L2_p = zeros(4,1); R_p = zeros(4,1); eoc_p=zeros(3,1); Reoc_p=zeros(3,1); conv=zeros(4,S); 
iit = [5 10 20 40 80 160 320 640]; 
ii=4;
for i=1:ii 
    it=iit(i);
    dx=1/it;
    str=['dx=',num2str(dx,'%.2e')];
    strvect(i,1:length(str))=str;
    [p,pE,q,qE,eps] = R2D_Non_steady(it,phis,E); % ,eps_err
    u=chebfun2(q,[0 1 0 1],'equi');
    uE=chebfun2(qE,[0 1 0 1],'equi');
    phisE=phi0(uE);      % source term as chebfun
    f=phi0(u)-phisE;     % modified eq. for man. sol.
    % f=phi0(u)-phis;     % modified eq. for man. sol.
    rezu=norm(f);
    L2_p(Itest) = ( dx )^(1/2) *norm(p-pE);
    R_p(Itest) = rezu; 
    if Itest >1
        eoc_p(Itest-1)=log10(L2_p(Itest-1)/L2_p(Itest))/log10(2);
        Reoc_p(Itest-1)=log10(R_p(Itest-1)/R_p(Itest))/log10(2);
    end
    Itest = Itest+1;
    conv(i,:)=eps;
end
figure;
for i=1:ii 
    P(i)=plot(iS(1:10:end),conv(i,1:10:end)); hold on;
end
mark={'o','+','x','*','d','^','v','s'}';
NameArray = {'Marker'}; ValueArray = mark(1:ii); 
set(P,NameArray,ValueArray);
set(gca,'yscale','log'); box on;
xlabel('$s$','Interpreter','latex');
ylabel('$\|\psi^s - \psi_{ref}\|$','Interpreter','latex');
legend(strvect,Location="best"); legend('boxoff');

disp('Error norms');
fprintf('||p-pE|| = %0.2e \n',L2_p)
disp('Order of Accuracy');
fprintf('EOA_err : %0.2e \n',eoc_p)

disp('Residuals');
fprintf('||phi(u)|| = %0.2e \n',R_p)
disp('Order of Accuracy');
fprintf('EOA_res = %0.2e \n',Reoc_p)

p1=pE;
figure;
xx=0:dx:1;
plot(xx,p1,'-*',xx,p,'-+');
xlabel('$z$','Interpreter','latex');
ylabel('$\psi(z,t)$','Interpreter','latex');
legend('analytical','numerical','Location','best'); legend('boxoff')

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,pE,q,qE,eps] =  R2D_Non_steady(it,phis,E)
%%   Grid Initialization
I = it+1;
a=0; b=1;
dx = (b-a)/(I-1);
x = a:dx:b;
%%   Parameters
T=1;
past=T/3;
S=10000;
L=5; % 1; % 
maxr=0.5;% 0.8; % 
Tolerance = 1e-6;
%%  Parametrization functions
theta = @(c)  1./(1-c);
K = @(c) c.^2;
%%   Initial Conditions
p0=E(0,:); c0=p0(x); % exact solution given by the chebfun 'E' evaluated at 't=0' and '(X,Y)'
c = c0;
%% solution
t=0; kt=1; tt=0;
eps=zeros(1,S); 
cc=zeros(1,I);
tht = theta(c);
tht0=tht;
ca=c;
while t<=T
DK=K(c);
D=(DK(1:I-1)+DK(2:I))/2;
dt=L*dx^2*maxr/max(D)/2;
t=t+dt; tt=tt+1;
F=phis(t,:); fs=F(x);
    for s=1:S

        %% FD scheme:
        DK=K(c);
        D=(DK(1:I-1)+DK(2:I))/2;
        r=D*maxr/(2*max(D)); 
        rloc=1-(r(1:I-2)+r(2:I-1));
        cc(2:I-1)=rloc.*c(2:I-1)+r(2:I-1).*c(3:I)+r(1:I-2).*c(1:I-2);
        %% boundary conditions
        %%%% BC_Left/Right
        cc(1)=c0(1);
        cc(I)=c0(I);
        %% Source term
        dtht=(tht0-tht)/L;
        f=diff(r)*dx+dtht(2:I-1)+fs(2:I-1)*dt/L;
        cc(2:I-1)=cc(2:I-1)+f;
        c=cc;
        %% Convergence criterion
        tol_eps=norm(c-ca); 
        eps(s)=tol_eps;
        if tol_eps <= Tolerance
            break
        end
        tht=theta(c);
        ca=c;
    end
    if  kt*past>=t && kt*past<t+dt && t<=T %tt==TF %
        rndt=round(kt*past,2);
        % fprintf('kt*past= %d\n',rndt);
        kt=kt+1;
    end
    tht0=tht;
    %%
    q(:,tt)=c; % numerical solution at time step 'tt' ===> full numerical solution in space-time domain
    e=E(t,:); qE(:,tt)=e(x)';
    % qE(:,tt)=E(t,x)'; 
    %%
end
e=E(t,:);
pE=e(x);         % evaluated chebfun
p=c;
fprintf('dx = %0.2e   dt = %0.2e\n',dx,dt);
end

%% Accuracy_Chebfun2
% dx = 2.00e-01   dt = 9.02e-02
% dx = 1.00e-01   dt = 2.23e-02
% dx = 5.00e-02   dt = 5.56e-03
% dx = 2.50e-02   dt = 1.39e-03
% Error norms
% ||p-pE|| = 1.50e-02 
% ||p-pE|| = 3.61e-03 
% ||p-pE|| = 8.61e-04 
% ||p-pE|| = 1.32e-04 
% Order of Accuracy
                                    % EOA_err : 2.05e+00 
                                    % EOA_err : 2.07e+00 
                                    % EOA_err : 2.70e+00 
% Residuals
% ||phi(u)|| = 8.78e-02 
% ||phi(u)|| = 1.93e-02 
% ||phi(u)|| = 4.57e-03 
% ||phi(u)|| = 9.73e-04 
% Order of Accuracy
                                    % EOA_res = 2.18e+00 
                                    % EOA_res = 2.08e+00 
                                    % EOA_res = 2.23e+00 
% Elapsed time is 12.322953 seconds.
