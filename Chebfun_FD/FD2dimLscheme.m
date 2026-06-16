%% Two dimensional FD L-scheme code to solv Richards equation

clear all; close all
tic
%%   Grid Initialization
I=41; J = I;
a=0; b=1;
c=0; d=1;
dx = (b-a)/(I-1);
x = a:dx:b;
dy=(d-c)/(J-1);
y=c:dy:d;
%%   Parameters
T=1;
past=T/3;
S=1000;
L=1; 
maxr=0.8;
Tolerance = 1e-5;
%%  Parametrization functions
theta = @(c)  1./(1-c);
K = @(c) c.^2; %
solE = @(t,x,y)  -t .* x .* (1-x) .* y .* (1-y)-1;
F =  @(t,x,y) -x.*(1-x).*y.*(1-y)./(2+t.*x.*(1-x).*y.*(1-y)).^2+...
    -2*t*(1+t.*x.*(1-x).*y.*(1-y)).^2.*(x.*(1-x)+y.*(1-y)) ...
    +2*(1+t.*x.*(1-x).*y.*(1-y)).*((-t*(1-2*x).*y.*(1-y)).^2-t*(-t*x.*(1-x).*(1-2*y)+1).*x.*(1-x).*(1-2*y));
%%   Initial Conditions
[X,Y] = meshgrid(x,y) ;
c0=solE(0,X,Y);
sumcinit=sum(sum(c0));
c = c0;
%% solution
tgraf=0; tconv=0; t=0; iT=0; kt=1;
sumc=zeros(1,floor(T/past)); sumn=zeros(1,floor(T/past)); tvect=zeros(1,floor(T/past));
eps=zeros(1,S); cc=zeros(J,I);
tevol=zeros(1,floor(T/past)); cgraf=zeros(floor(T/past),J,I); thtgraf=zeros(floor(T/past),J,I);
tht = theta(c);
tht0=tht;
ca=c;
iS=1:S;
DK=K(c);
Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
D=DK(2:J-1,2:I-1);
dt=2*(max(max(Dx))/dx^2/L+max(max(Dy))/dy^2/L); dt=maxr*1/dt; 
t=dt;
while t<=T
    fs=F(t,X,Y);
    iT=iT+1;
    for s=1:S
        DK=K(c);
        Dx=(DK(2:J-1,1:I-1)+DK(2:J-1,2:I))/2;
        Dy=(DK(1:J-1,2:I-1)+DK(2:J,2:I-1))/2;
        D=DK(2:J-1,2:I-1);
        rx=dt*Dx/dx^2/L; ry=dt*Dy/dy^2/L; 
        rloc=1-(rx(:,1:I-2)+rx(:,2:I-1)+ry(1:J-2,:)+ry(2:J-1,:));
        cc(2:J-1,2:I-1)=rloc.*c(2:J-1,2:I-1) ...
            +rx(:,1:I-2).*c(2:J-1,1:I-2)+rx(:,2:I-1).*c(2:J-1,3:I) ...
            +ry(1:J-2,:).*c(1:J-2,2:I-1) +ry(2:J-1,:).*c(3:J,2:I-1);
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
        if kt*past>=t && kt*past<t+dt && t<=T %tt==TF %
            eps(s)=tol_eps;
        end
        if tol_eps <= Tolerance
            break
        end
        tht=theta(c);
        ca=c;
    end
    if  kt*past>=t && kt*past<t+dt && t<=T %tt==TF %
        tgraf=tgraf+1;
        rndt=round(kt*past,2);
        str=['t=',num2str(rndt)];
        strvect(tgraf,1:length(str))=str;
        figure(1); box; hold on;
        P(tgraf)=plot(iS,eps);
        tvect(tgraf)=rndt;
        kt=kt+1;
    end
    tht0=tht;
    tE=t;
    t=t+dt;
    %%
    tt=round((t-dt)/dt);
    q(tt,:,:)=c;
end
sumcfin=sum(sum(c));
%% Velocity components
Vx=-D.*((c(2:J-1,3:I)-c(2:J-1,1:I-2))/(2*dx));
Vy=-D.*((c(3:J,2:I-1)-c(1:J-2,2:I-1))/(2*dy)+1);
%% plots
NameArray = {'Marker'}; ValueArray = {'o','+','x'}';
set(P,NameArray,ValueArray);
set(gca,'yscale','log'); box on;
xlabel('$s$','Interpreter','latex');
ylabel('$\|\psi^s - \psi^{s-1}\|$','Interpreter','latex');
legend(strvect); legend('boxoff');
figure(2);
mesh(x,y,c);
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
zlabel('$\psi(x,z,t)$','Interpreter','latex'); view(115,15); %(20,50);
grid on ;
cE=solE(tE,X,Y);
figure(3);
mesh(x,y,c-cE);
fprintf('norm(c-cE) =  %0.2e \n',norm(c-cE)) ;
xlabel('$x$','Interpreter','latex'); ylabel('$z$','Interpreter','latex');
zlabel('$\psi(x,z,t)-\psi(x,z,t)_{exact}$','Interpreter','latex'); view(115,15); %(20,50);
grid on ;
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
fprintf('The final time is : %0.2e \n',t-dt) ;
toc
% Elapsed time is 5.073084 seconds (I = 11;)
% Elapsed time is 70.044316 seconds (for I = 41;)