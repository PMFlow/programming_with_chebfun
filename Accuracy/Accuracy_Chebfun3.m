%%   Estimated order of accuracy for NON-steady Richards-2D

close all; clear all;
tic
T=1;
x1=0; x2=1;
y1=0; y2=1;
D=[0 T x1 x2 y1 y2];
x = chebfun3(@(x,y,z) x,D);
y = chebfun3(@(x,y,z) y,D);
z = chebfun3(@(x,y,z) z,D);
E=chebfun3(@(x,y,z) -x.*y.*(1-y).*z.*(1-z)-0.5,[0 1 0 1 0 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi0 = @(r) diffx(r)./((1-r).^2)-r.^2.*(diffy(r,2)+diffz(r,2))-2*r.*(diffy(r).*diffy(r)+diffz(r).*(diffz(r)+1));
phis=phi0(E);      % source term as chebfun
Itest = 1;
S= 10000; iS=1:S;
L2_p = zeros(4,1); R_p = zeros(4,1); eoc_p=zeros(3,1); Reoc_p=zeros(3,1); conv=zeros(4,S); rezu=0; %conv_err=zeros(4,S);
iit = [5 10 20 40 80 160 320 640]; 
ii=3;%4;
for i=1:ii
    it=iit(i);
    dx=1/iit(i); dy=dx;
    str=['dx=',num2str(dx,'%.2e')];
    strvect(i,1:length(str))=str;
    [p,pE,q,eps] = R3D_Non_steady(it,phis,E); % ,eps_err
    u3=chebfun3(q,[0 T 0 1 0 1],'vectorize','equi'); % ,'vectorize','splitting on'
    f=phi0(u3)-phis;     % modified eq. for man. sol.
% %
%     uE=chebfun3(qE,[0 T 0 1 0 1],'equi');
%     phisE=phi0(uE);      % source term as chebfun
%     f=phi0(u3)-phisE;
% %
    rezu=norm(f);
    L2_p(Itest) = (dx*dy)^(1/2)*norm(p-pE);
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
    P(i)=plot(iS(1:100:end),conv(i,1:100:end)); hold on;
end
mark={'o','+','x','*','d','^','v','s'}';
NameArray = {'Marker'}; ValueArray = mark(1:i);
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

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,pE,q,eps] = R3D_Non_steady(it,phis,E)
%%   Grid Initialization
I = it+1; J = I;
a=0; b=1;
c=0; d=1;
dx = (b-a)/(I-1);
x = a:dx:b;
dy=(d-c)/(J-1);
y=c:dy:d;
%%   Parameters
T=1;
past=T/3;
S=10000;
% L=1200; % w. "-1"
% maxr=0.5;

% L=1000; maxr= 0.1; % w. "-0.5" - pt. colocviul F-R
% L=1800; maxr=0.1; % 0.085;
L=3000; maxr= 0.05; 
Tolerance = 1e-6;
%%  Parametrization functions
theta = @(c)  1./(1-c);
K = @(c) c.^2;
%%   Initial Conditions
[X,Y] = meshgrid(x,y) ;
p0=E(0,:,:); c0=p0(X,Y); % exact solution given by the chebfun3 'E' evaluated at 't=0' and '(X,Y)'
c = c0;
% figure; surf(p0); figure; surf(c0);
%% solution
t=0; kt=1; tt=0;
eps=zeros(1,S); cc=zeros(J,I);
tht = theta(c);
tht0=tht;
ca=c;
    DK=K(c);
    Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
    Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
    dt=12*(max(max(Dx))/(L*dx^2)+max(max(Dy))/(L*dy^2)); dt=maxr*1/dt;
    t=t+dt; tt=tt+1;
while t<=T
    DK=K(c);
    Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
    Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
    dt=12*(max(max(Dx))/(L*dx^2)+max(max(Dy))/(L*dy^2)); dt=maxr*1/dt;
    % t=t+dt; tt=tt+1;
    F=phis(t,:,:); fs=F(X,Y);
    for s=1:S
        DK=K(c);
        Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
        Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
        D=DK(2:J-1,2:I-1);
        rx=dt*Dx/dx^2/L; ry=dt*Dy/dy^2/L;
        rloc=1-(rx(:,1:I-2)+rx(:,2:I-1)+ry(1:J-2,:)+ry(2:J-1,:));
        cc(2:J-1,2:I-1)=rloc.*c(2:J-1,2:I-1) ...
            +rx(:,1:I-2).*c(2:J-1,1:I-2)+rx(:,2:I-1).*c(2:J-1,3:I) ...
            +ry(1:J-2,:).*c(1:J-2,2:I-1)+ry(2:J-1,:).*c(3:J,2:I-1);
        %% boundary conditions
        %%%% BCXLeft/Right
        cc(:,1)=c0(:,1);
        cc(:,I)=c0(:,I);
        %%%% BCYBottom/Upper
        cc(1,:)=c0(1,:);
        cc(J,:)=c0(J,:);
        %% Source term
        dtht=(tht0-tht)/L;
        f=(ry(2:J-1,:)-ry(1:J-2,:))*dy + dtht(2:J-1,2:I-1)...
            + fs(2:J-1,2:I-1)*dt/L;
        cc(2:J-1,2:I-1)=cc(2:J-1,2:I-1)+f;
        c=cc;
        %% Convergence criterion
        tol_eps=norm(c-ca);
        if kt*past>=t && kt*past<t+dt && t<=T
            eps(s)=tol_eps;
        end
        if tol_eps <= Tolerance
            break
        end
        ca=c;  tht=theta(c);
    end
    if  kt*past>=t && kt*past<t+dt && t<=T %tt==TF %
        rndt=round(kt*past,2);
        % fprintf('kt*past = %d\n',rndt);
        kt=kt+1;
    end
    % tht=theta(c);
    tht0=tht;
    %%
    q(tt,:,:)=c; % numerical solution at time step 'tt' ===> full numerical solution in space-time domain
    % % e=E(tt,:,:); 
    % % qE(tt,:,:)=e(X,Y); 
    % e=E(t,:,:); 
    % qE(:,:,tt)=e(X,Y); 
    %%
    % fprintf('t = %d\n',t);
    t=t+dt; tt=tt+1;
end
e=E(t,:,:);
pE=e(X,Y);         % evaluated chebfun
p=c;
fprintf('dx = %0.2e   dt = %02e\n',dx,dt);
end

%% corrected
%% Accuracy_Chebfun3: L=3000; maxr= 0.05; 
% dx = 2.00e-01   dt = 1.000000e+00
% dx = 1.00e-01   dt = 2.107507e-01
% dx = 5.00e-02   dt = 5.052996e-02
% dx = 2.50e-02   dt = 1.247695e-02
% Warning: Matrix is close to singular or badly scaled. Results may be
% inaccurate. RCOND =  1.572017e-16. 
% Warning: Matrix is close to singular or badly scaled. Results may be
% inaccurate. RCOND =  1.801698e-16. 
% Error norms
% ||p-pE|| = 3.30e-02 
% ||p-pE|| = 7.27e-03 
% ||p-pE|| = 2.13e-03 
% ||p-pE|| = 1.35e-03 
% Order of Accuracy
                                    % EOA_err : 2.18e+00 
                                    % EOA_err : 1.77e+00 
                                    % EOA_err : 6.55e-01 
% Residuals
% ||phi(u)|| = 1.17e-01 
% ||phi(u)|| = 2.61e-02 
% ||phi(u)|| = 7.70e-03 
% ||phi(u)|| = 4.72e-03 
% Order of Accuracy
                                    % EOA_res = 2.16e+00 
                                    % EOA_res = 1.76e+00 
                                    % EOA_res = 7.06e-01 
% Elapsed time is 147.072107 seconds.


