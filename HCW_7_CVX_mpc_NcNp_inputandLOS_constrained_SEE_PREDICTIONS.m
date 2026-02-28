% Clohessy Wiltshere relative dynamics model 
% Difference between the previous one is, in this code we are able to see
% the predictions of MPC in our plots.
% Proximity operation with MPC via CVX
% "Reverse Engineering" by ISKENDER.O.B, 02/27/2017 
% Supervisor: Assoc. Prof. LING Keck Voon
clear ;
close all;
clc;
%==========================================================================
% Variable Name         Description           Units
%==========================================================================
%      X0        Initial Condition Vector      [m m m m/s m/s m/s]
%      x0        Initial Relative X Positon    km
%      y0        Initial Relative Y Position   km
%      z0        Initial Relative Z Position   km
%      x0dot     Initial Relative X Velocity   km/s
%      y0dot     Initial Relative Y Velocity   km/s
%      z0dot     Initial Relative Z Velocity   km/s
%      mu        Gravitational parameter       m^3/s^2
%      a         Semi Major Axis               m
%      tf        Time Final Position Time      min
%      t         Simulation Time               sec
%      Ts        Sampling Time                 sec
%      n         Mean Motion                   Rad/sec
%      Q         State Weight Matrix
%      R         Input Weight Matrix
%      Qr        Terminal penalty Weight Matrix
%      Np         Prediction Horizon
%      Nc        Control Horizon
%      Q         State Weight Matrix
%      R         Input Weight Matrix
%      Qp        Terminal Penalty Weight Matrix
%==========================================================================
%Assign Model and Controller parameters;
mpc.t        = 200;
% mpc.X0       = [1500 1000 1000 0 0  0]'; 
% mpc.X0       = [0 100 0 1 1 -1]'; 
mpc.X0       = [0 200 0 4 0 -2]'; 

mpc.Np       = 50;  
mpc.Nc       = 50;  
mpc.Ts       = 1;  
mpc.mu       = 398600;
mpc.a        = 6968.1363;
mpc.omega    = sqrt(mpc.mu/(mpc.a^3));
mpc.States_n = 6;
mpc.Inputs_m = 3;
mpc.u_max    =  0.1;
mpc.u_min    = -0.1;
mpc.Q        = eye(6);
% mpc.R        = diag([2*1e8 1.6*1e8 1e8]);
% mpc.R        = diag([2*1e7 1.6*1e7 1e7]);
% mpc.R        = diag([2 1.6 1]);
mpc.R = 1e4*eye(3);
mpc.A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     3*mpc.omega^2 0 0 0 2*mpc.omega 0;
     0 0 0 -2*mpc.omega 0 0;
     0 0 -mpc.omega^2 0 0 0];
mpc.B = [zeros(3); eye(3)];
mpc.C = eye(6);
mpc.D = zeros(6,3);
%Create CTS time State Space Model
mpc.ctsModel = ss(mpc.A,mpc.B,mpc.C,mpc.D);
%Transform CTS time State Space Model to Discrete Model
mpc.dModel   = c2d(mpc.ctsModel,mpc.Ts);

[mpc.K,mpc.Qp]   = lqry(mpc.dModel,mpc.Q,mpc.R);
[mpc.Ad, mpc.Bd] = c2dm(mpc.A,mpc.B,mpc.C,mpc.D,mpc.Ts);

%Initialize the data to be saved
mpc.u_final = [];
mpc.XX      = mpc.X0;
mpc.Tcompt  = [];

%Docking Port Dimention & Slope of the LOS constraint Lines
mpc.LOS.x_0 = 1.5 ;% docking port dimention
mpc.LOS.z_0 = 1.5 ;% docking port dimention
mpc.LOS.c_x = 1 ;  %slope of the tetrahedral region
mpc.LOS.c_z = 1 ;  %slope of the tetrahedral region

mpc.LOS.A_c =[0 -1 0 0 0 0 ; ...
                   mpc.LOS.c_x -1 0 0 0 0; ...
                  -mpc.LOS.c_x -1 0 0 0 0; ...
                   0 -1 mpc.LOS.c_z 0 0 0; ...
                   0 -1 -mpc.LOS.c_z 0 0 0];

mpc.LOS.b_c = [0 mpc.LOS.c_x*mpc.LOS.x_0 mpc.LOS.c_x*mpc.LOS.x_0 mpc.LOS.c_z*mpc.LOS.z_0 mpc.LOS.c_z*mpc.LOS.z_0]';
mpc.XX0 = mpc.X0;  
mpc.OOO = [];
mpc.KKK = [];
mpc.MMM = [];
for i = 1: mpc.t
    ix = 0;
        fprintf('\n \n \n \n                   Elapsed time  = %d [s] \n \n \n',i)  %Knowing that each step is 1s
        Distance = sqrt(mpc.X0(1)^2+mpc.X0(2)^2+mpc.X0(3)^2);
        fprintf('                   Distance from the target = %.4f [metre]\n ',Distance)  %Knowing that each step is 1s
            mpc.myJ = 0;
       tic
 cvx_begin quiet
      variables X(mpc.States_n,mpc.Np+1) U(mpc.Inputs_m,mpc.Np)
if mpc.Np < mpc.Nc
    fprintf('\n\n\nControl Horizon(mpc.Nc) can not be greater than Prediction Horizon(mpc.Np)!!!! \n\nLook at your mistake: Np %d <?? Nc %d \n',mpc.Nc,mpc.Np)
   return
else
    for ix = 1:(mpc.Np)
        if(ix<=mpc.Nc)        
            mpc.myJ = mpc.myJ +  (X(:,ix))'*mpc.Q*(X(:,ix)) + (U(:,ix))'*mpc.R*(U(:,ix));
        else   % If Nc is less than Np
            U(:,ix) == U(:,mpc.Nc);
            mpc.myJ = mpc.myJ +  (X(:,ix))'*mpc.Q*(X(:,ix)) + (U(:,ix))'*mpc.R*(U(:,ix));
        end
    end
            mpc.myJ = mpc.myJ + (X(:,mpc.Np+1))'*mpc.Qp*(X(:,mpc.Np+1));
 
end
        X(:,2:mpc.Np+1)     == mpc.Ad*X(:,1:mpc.Np) + mpc.Bd*U;
        X(:,1)              == mpc.X0;       
         mpc.z0            =  X(:,:);
        %Input Constraint Formulation
   - mpc.u_max <= U(:,:) <= mpc.u_max
hold on
        %Line of Sight Constraint Formulation
 for ii= 1:mpc.Np+1
    mpc.LOS.A_c*X(:,ii)  <= mpc.LOS.b_c;
 end

            cvx_solver SDPT3
%             cvx_solver Mosek
%             cvx_solver Sedumi
%             cvx_solver Gurobi

            minimize mpc.myJ
    cvx_end   
    El_time = toc;
    mpc.Tcompt = [mpc.Tcompt; El_time];   
    mpc.u = U(:,1);
    mpc.xxx = X(1,:);
    mpc.yyy = X(2,:);
    mpc.zzz = X(3,:);   
    mpc.OOO = [mpc.OOO; mpc.xxx];
    mpc.KKK = [mpc.KKK; mpc.yyy];
    mpc.MMM = [mpc.MMM; mpc.zzz];
    mpc.X0 = mpc.Ad*mpc.X0+mpc.Bd*mpc.u;
    mpc.XX = [mpc.XX mpc.X0]; 
    mpc.u
    mpc.X0
    mpc.u_final = [mpc.u_final mpc.u];
end
fprintf('End of the simulation\n')
u_final = mpc.u_final;
XX      = mpc.XX;
Tcompt  = mpc.Tcompt;

figure(1)
subplot(3,1,1)
plot(XX(1,:))
set(gca,'XLim',[1 mpc.t+1])
grid on
subplot(3,1,2)
plot(XX(2,:))
set(gca,'XLim',[1 mpc.t+1])
grid on
subplot(3,1,3)
plot(XX(3,:))
set(gca,'XLim',[1 mpc.t+1])
grid on

figure(2)
subplot(3,1,1)
plot(u_final(1,:))
set(gca,'XLim',[1 mpc.t+1])
grid on
subplot(3,1,2)
plot(u_final(2,:))
set(gca,'XLim',[1 mpc.t+1])
grid on
subplot(3,1,3)
plot(u_final(3,:))
set(gca,'XLim',[1 mpc.t+1])
grid on
%%

figure(4)
hold on
plot3(mpc.XX0(1),mpc.XX0(2),mpc.XX0(3),'y.','MarkerSize',36)

h1 = plot3(XX(1,:),XX(2,:),XX(3,:),'r*','MarkerSize',3);
axis tight
% axis([-0.2 0.5 -5 60 -1 1 ])
% axis([-100 100 -5 1100 -1 1 ])

%%
if exist('XX1') ==1
h2 = plot3(XX1(1,:),XX1(2,:),XX1(3,:),'bd','MarkerSize',5);
end
% axis([-1000 1000 -1100 1100 -1000 1000 ])

% axis([-100 100 -5 1100 -1 1 ])


%%
hold on
plot3(0,0,0,'g.','MarkerSize',36)
% plot3(0,-250,100,'y.','MarkerSize',36)
L = 100;
axis tight
axis([-100 100 -1100 110 -1000 1000 ])
axis([-300 300 -300 300 -300 1000 ])

c = 1;

x0 = 1.5;
z0 = 1.5;

x  = 1:L;
y  = x-1.5;
z  = x;
h3 = plot3(x,y,z,'k.','MarkerSize',11);
xl = xlabel('Radial, x');
yl = ylabel('In-track, y');
zl = zlabel('Cross-track, z');
hold on
x1 = (1:L);
y1 = -x-x0;
z1 = x1;
% plot3(x1,y1,z1)

x1 = -1:-1:-L;
y1 = -x1-x0;
z1 = x1;
h4 = plot3(x1,y1,z1,'k.','MarkerSize',11);

x2 = 1:L;
y2 = x2-x0;
z2 = -x2;
h5 = plot3(x2,y2,z2,'k.','MarkerSize',11);
% 
% 
x3 = -1:-1:-L;
y3 = -x3-x0;
z3 = -x3;
plot3(x3,y3,z3,'k.','MarkerSize',11);
% plot3(mpcOBJ.X0(1),mpcOBJ.X0(2),mpcOBJ.X0(3),'k.','MarkerSize',11)

axis([-200 200 -110 1100 -1000 1000 ]);
grid on
% az = -58; %Azimuth
% el = 17; %Elevation
az = 0; %Azimuth
el = 90; %Elevation
% az = -90 %Azimuth
% el = 0 %Elevation
view(az, el);
axis tight

title1 = title(['Line of Sight Constraint Comparison']);
grid on
leg1 = legend('Without LOS','With LOS','Docking Region','Initial Position','Location','Best');
% leg1 = legend('Without LOS','Location','Docking Region','Initial Position','Best');
% leg1 = legend('Without LOS','Docking Region','Initial Position','Location','Best');
% leg1 = legend('Docking Region','Initial Position','Location','Best');

set([title1 xl yl zl leg1],'interpreter','latex','fontsize',22)
set([leg1],'interpreter','latex','fontsize',16)
axis([-1000 1000 -1000 1000 -1000 1000 ])

%%
figure(61)
OOO = mpc.OOO;
KKK = mpc.KKK;
MMM = mpc.MMM;
for ii = 1:mpc.t
plot(ii:ii+mpc.Np,OOO(ii,:))
hold on
end
hold on
plot(mpc.XX(1,:),'r.','MarkerSize',14)
figure(62)
for ii = 1:mpc.t
plot(ii:ii+mpc.Np,KKK(ii,:))
hold on
end
hold on
plot(mpc.XX(2,:),'r.','MarkerSize',14)
figure(63)
for ii = 1:mpc.t
plot(ii:ii+mpc.Np,MMM(ii,:))
hold on
end
hold on
plot(mpc.XX(3,:),'r.','MarkerSize',14)
%%
figure (3)
% plot3(mpc.XX0(1),mpc.XX0(2),mpc.XX0(3),'y.','MarkerSize',36)
for ii = 1:mpc.t
% plot3(OOO(ii,:),KKK(ii,:),MMM(ii,:),'b--','MarkerSize',16)
% plot3(OOO(ii,:),KKK(ii,:),MMM(ii,:),'k--','MarkerSize',16)
plot3(OOO(ii,:),KKK(ii,:),MMM(ii,:),'r--','MarkerSize',16)

hold on
end




% figure(4)
hold on

% h1 = plot3(XX(1,:),XX(2,:),XX(3,:),'r*','MarkerSize',7);
h1 = plot3(XX(1,:),XX(2,:),XX(3,:),'yd','MarkerSize',7);
h1 = plot3(XX(1,:),XX(2,:),XX(3,:),'md','MarkerSize',7);

axis tight
% axis([-0.2 0.5 -5 60 -1 1 ])
% axis([-100 100 -5 1100 -1 1 ])

%%
if exist('XX1') ==1
h2 = plot3(XX1(1,:),XX1(2,:),XX1(3,:),'bd','MarkerSize',5);
end
% axis([-1000 1000 -1100 1100 -1000 1000 ])

% axis([-100 100 -5 1100 -1 1 ])


%%
hold on
plot3(0,0,0,'g.','MarkerSize',36)
% plot3(0,-250,100,'y.','MarkerSize',36)
L = 100;
axis tight
axis([-100 100 -1100 110 -1000 1000 ])
axis([-300 300 -300 300 -300 1000 ])

c = 1;

x0 = 1.5;
z0 = 1.5;

x  = 1:L;
y  = x-1.5;
z  = x;
h3 = plot3(x,y,z,'k.','MarkerSize',11);
xl = xlabel('Radial, x');
yl = ylabel('In-track, y');
zl = zlabel('Cross-track, z');
hold on
x1 = (1:L);
y1 = -x-x0;
z1 = x1;
% plot3(x1,y1,z1)

x1 = -1:-1:-L;
y1 = -x1-x0;
z1 = x1;
h4 = plot3(x1,y1,z1,'k.','MarkerSize',11);

x2 = 1:L;
y2 = x2-x0;
z2 = -x2;
h5 = plot3(x2,y2,z2,'k.','MarkerSize',11);
% 
% 
x3 = -1:-1:-L;
y3 = -x3-x0;
z3 = -x3;
plot3(x3,y3,z3,'k.','MarkerSize',11);
% plot3(mpcOBJ.X0(1),mpcOBJ.X0(2),mpcOBJ.X0(3),'k.','MarkerSize',11)

axis([-200 200 -110 1100 -1000 1000 ]);
grid on
% az = -58; %Azimuth
% el = 17; %Elevation
az = 0; %Azimuth
el = 90; %Elevation
% az = -90 %Azimuth
% el = 0 %Elevation
view(az, el);
axis tight

title1 = title(['Line of Sight Constraint Comparison']);
grid on
leg1 = legend('Real State','Reference','Predicted States','Location','Best');

% leg1 = legend('Without LOS','With LOS','Docking Region','Initial Position','Location','Best');
% leg1 = legend('Without LOS','Location','Docking Region','Initial Position','Best');
% leg1 = legend('Without LOS','Docking Region','Initial Position','Location','Best');
% leg1 = legend('Docking Region','Initial Position','Location','Best');

set([title1 xl yl zl leg1],'interpreter','latex','fontsize',22)
set([leg1],'interpreter','latex','fontsize',16)
axis([-1000 1000 -100 1000 -1000 1000 ])
%%
save('C:\Users\ISKE0001\Dropbox\LOS_NOLOS\Sava_data\Constrained_vs_unconstrainedMPC\sim1_Thales_sunum.mat','XX','mpc')