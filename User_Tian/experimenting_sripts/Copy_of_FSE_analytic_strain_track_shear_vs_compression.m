%Copy_of_FSE_analytic_strain_track_shear_vs_compression.m
% Finite Strain Ellipse from different combination of strain history
%clear
clc
close all
clear

%% asign location and size of figures
%ref: https://www.mathworks.com/matlabcentral/answers/684798-how-to-display-figure-on-a-specific-monitor-when-there-are-multiple-screens
%ref: https://www.mathworks.com/help/releases/R2018b/matlab/ref/matlab.ui.root-properties.html#buc8_0n-MonitorPositions
p = get(0, "MonitorPositions"); 
%p(1,1) is origin X0
%p(1,2) is origin Y0
%p(1,3) is width of primary display
%p(1,4) is height of primary display
%
%p(2,1) is origin p2X0 for external display relative to X0
%p(2,2) is origin p2Y0 for external display relative to Y0
%p(2,3) is width of external display
%p(2,4) is height of external display

p1X0=p(1,1); %is origin X0
p1Y0=p(1,2); % is origin Y0
p1W=p(1,3);  % is width of primary display p1
p1H=p(1,4);  % is height of primary display p2

p2X0=p(2,1); % is origin p2X0
p2Y0=p(2,2); % is origin p2Y0
p2W=p(2,3);  % is width of external display p2
p2H=p(2,4);  % is height of external display p2

f1 = figure(1);
p2_scale = .3;
f1_position = [p2X0,p2Y0,p2_scale*p2W,p2_scale*p2H];
f1.Position = f1_position;


%% Initialization
epsilon_dot = 1e-15;% 1/s
dvxdy = epsilon_dot;
%% Imposed velocity gradient
%L = [0 dvxdy; 0 0] % simple shear
L = [0 0; dvxdy 0] % simple shear
%L = [0 dvxdy; -dvxdy 0] % rigid rotation counterclockwise 
%L = [0 dvxdy; dvxdy 0] % L pure shear
%L = [dvxdy 0; 0 -dvxdy] % L compression vertically
%L = [dvxdy dvxdy; -dvxdy -dvxdy] % L compression vertically + simple shear

time_end = 50; %100; % Myrs
sec_in_yr = 365.25*24*3600; %seconds
Ntime = 200;
%Ntime = 500;
t = linspace(0,time_end*1e6*sec_in_yr,Ntime);
dt = t(2)-t(1);
finite_strain = zeros(1,Ntime);
angle_from_x = zeros(1,Ntime);

F = [1 0; 0 1]; % initial deformation gradient tensor
%F = [1 0; 0 1.5]; % initial deformation gradient tensor
%F = [2 0; -.1 1];
%F = [2.2 0; 0 1]; % initial deformation gradient tensor

%F = [2 0; 0 3]; % initial deformation gradient tensor
%F = [1 0; 0 3]; % initial deformation gradient tensor

% getting FSE info from the initial F(0)
[FSEa,FSEb,FSEeigV,FSEtheta,FSEQ]=FSE(F);
FSEscale = 1; % make them look bigger

%% plotting the FSE with init F(time=0)
figure(f1)
ellipse(FSEscale*FSEa,FSEscale*FSEb,FSEtheta,0, ...
0,[0.9290 0.6940 0.1250],100,12); hold on;

axis equal;
grid on;
xlabel("X");
ylabel("Y");
set(gca, 'YDir','reverse');
fontsize(22,"points");

%% set up the Lq for tracking sense of deformation relative to FSE
Lq = zeros(2,2,Ntime); % Lq = QLQ'; strain rate relative to a-b coords.
finite_strain_ab = zeros(1,Ntime); % finite strain from FSE(a-b), directly QFQ'
angle_from_a = zeros(1,Ntime);
finite_strain_q = zeros(1,Ntime); % finite strain from FSE(a-b) due to Lq
angle_from_a_q = zeros(1,Ntime);
%% pure shear Lqps
Lqps = zeros(2,2,Ntime); %pure shear component of Lq
finite_strain_qps = zeros(1,Ntime);
angle_from_a_qps = zeros(1,Ntime);
%% simple shear Lqss
Lqss = zeros(2,2,Ntime); %simple shear component of Lq
finite_strain_qss = zeros(1,Ntime);
angle_from_a_qss = zeros(1,Ntime);
%% rigid rotation Lqrr
Lqrr = zeros(2,2,Ntime); %rigid rotation component of Lq
finite_strain_qrr = zeros(1,Ntime);
angle_from_a_qrr = zeros(1,Ntime);

%% for loop to update F(t) and corresponding FSE(t)
for i=1:Ntime
    %i
    % F = F+dt*L*F; % update F in the end.
    [FSEa,FSEb,FSEeigV,FSEtheta,FSEQ]=FSE(F);
    
    %FSEtheta/pi*180;
    FSEscale = 1; % make FSE them look bigger
    if (rem(i,10)==1 )
        ellipse(FSEscale*FSEa,FSEscale*FSEb,FSEtheta,0,0,'k',100,1); hold on;
        quiver(0,0,FSEa*FSEeigV(1,2),FSEa*FSEeigV(2,2),'off','red'); % the 'off' before 'red' is for turning off auto scaling
        quiver(0,0,FSEb*FSEeigV(1,1),FSEb*FSEeigV(2,1),'off','blue'); 
    end
    if (i==Ntime )
        ellipse(FSEscale*FSEa,FSEscale*FSEb,FSEtheta,0,0,'k',100,10); hold on;
        quiver(0,0,FSEa*FSEeigV(1,2),FSEa*FSEeigV(2,2),'off','red',LineWidth=10); % the 'off' before 'red' is for turning off auto scaling
        quiver(0,0,FSEb*FSEeigV(1,1),FSEb*FSEeigV(2,1),'off','blue',LineWidth=10); 
    end
    finite_strain(i) = log(FSEa/FSEb);
    angle_from_x(i) = FSEtheta/pi*180;
    %ref:https://math.hmc.edu/funfacts/area-of-an-ellipse/
    area_FSE = pi*FSEa*sqrt(FSEeigV(1,2)^2+FSEeigV(2,2)^2) * FSEb*sqrt(FSEeigV(1,1)^2+FSEeigV(2,1)^2);
    %FSEa
    %FSEb
    title_text = sprintf("Finite strain=%.3f;    long axis angle=%.2f; area=%.3f (m^2) \n time=%.0f (kyrs);",...
        finite_strain(i),angle_from_x(i),area_FSE,t(i)/sec_in_yr/1e3);
    title(title_text);
    
    %% for plotting the FSE relative to a-b coords in green, so it doesn't rotate at all.
    Fp = FSEQ*F*FSEQ';
    %Fp = FSEQ'*F*FSEQ;

    [FpFSEa,FpFSEb,FpFSEeigV,FpFSEtheta,FpQ]=FSE(Fp);
    if (rem(i,10)==1 || i==Ntime)
        ellipse(FSEscale*FpFSEa,FSEscale*FpFSEb,FpFSEtheta,0,0,'green',100,1); hold on;
    end
    % calculate finite strain from FSE(a-b)
    finite_strain_ab(i) = log(FpFSEa/FpFSEb);
    angle_from_a(i) = FpFSEtheta/pi*180;
    if (i==1 )
        fprintf("Initial Finite strain for green a-b FSE is %.3f\n",log(FpFSEa/FpFSEb));
    end
    if ( i==Ntime)
        fprintf("Eventual Finite strain for green a-b FSE is %.3f\n",log(FpFSEa/FpFSEb));
    end

  
    %% vel gradient Lq w.r.t a-b coords    (F(i+1) = F(i)+dt*L(i+1)*F(i))
    Lq(:,:,i) = FSEQ*L*FSEQ'; % rotate coordinate system from X-Y to a-b for strain rate/vel grad tensor to get Lq in a-b sys
    %Lq(:,:,i) = FSEQ'*L*FSEQ;
    Lqps(1,1,i) = Lq(1,1,i);
    Lqps(2,2,i) = Lq(2,2,i);
    Lqss(1,2,i) = Lq(1,2,i)+Lq(2,1,i);
    Lqrr(1,2,i) = -Lq(2,1,i);
    Lqrr(2,1,i) = Lq(2,1,i);

    %% update Fp relative to a-b coords
    %Fpq = Fp+dt*Lq(:,:,i)*Fp;
    Ltensor_pq = (2*eye(2)+dt*Lq(:,:,i))/(2*eye(2)-dt*Lq(:,:,i));
    Fpq = (Ltensor_pq)*Fp;
    [FpqFSEa,FpqFSEb,FpqFSEeigV,FpqFSEtheta,FpqQ]=FSE(Fpq);
    finite_strain_q(i) = log(FpqFSEa/FpqFSEb)-finite_strain(i);
    angle_from_a_q(i) = FpqFSEtheta/pi*180;


    %Fpqps = Fp+dt*Lqps(:,:,i)*Fp;
    Ltensor_pqps = (2*eye(2)+dt*Lqps(:,:,i))/(2*eye(2)-dt*Lqps(:,:,i));
    Fpqps = (Ltensor_pqps)*Fp;
    [FpqpsFSEa,FpqpsFSEb,FpqpsFSEeigV,FpqpsFSEtheta,FpqpsQ]=FSE(Fpqps);
    finite_strain_qps(i) = log(FpqpsFSEa/FpqpsFSEb)-finite_strain(i);
    angle_from_a_qps(i) = FpqpsFSEtheta/pi*180;


    %Fpqss = Fp+dt*Lqss(:,:,i)*Fp;
    Ltensor_pqss = (2*eye(2)+dt*Lqss(:,:,i))/(2*eye(2)-dt*Lqss(:,:,i));
    Fpqss = (Ltensor_pqss)*Fp;
    [FpqssFSEa,FpqssFSEb,FpqssFSEeigV,FpqssFSEtheta,FpqssQ]=FSE(Fpqss);
    finite_strain_qss(i) = log(FpqssFSEa/FpqssFSEb)-finite_strain(i);
    angle_from_a_qss(i) = FpqssFSEtheta/pi*180;

    %Fpqrr = Fp+dt*Lqrr(:,:,i)*Fp;
    Ltensor_pqrr = (2*eye(2)+dt*Lqrr(:,:,i))/(2*eye(2)-dt*Lqrr(:,:,i));
    Fpqrr = (Ltensor_pqrr)*Fp;
    [FpqrrFSEa,FpqrrFSEb,FpqrrFSEeigV,FpqrrFSEtheta,FpqrrQ]=FSE(Fpqrr);
    finite_strain_qrr(i) = log(FpqrrFSEa/FpqrrFSEb)-finite_strain(i);
    angle_from_a_qrr(i) = FpqrrFSEtheta/pi*180;


    %----------------------------------------
    pause(.01);
    %F = F+dt*L*F;
    Ltensor = (2*eye(2)+dt*L)/(2*eye(2)-dt*L);
    F = (Ltensor)*F;
end
hold off;

f2 = figure(2);
f2_position = [p2X0+p2_scale*p2W,p2Y0,p2_scale*p2W,p2_scale*p2H];
f2.Position = f2_position;
%figure(2)
plot(t/sec_in_yr/1e6,finite_strain,'k',linewidth=20); hold on;
plot(t/sec_in_yr/1e6,finite_strain_ab,'r',linewidth=10); hold on;
plot(t/sec_in_yr/1e6,finite_strain_q,'blue',linewidth=5); hold on;
legend("X-Y",'a-b','Lq increments');
xlabel("time (Myr)")
ylabel("finite strain log(a/b)")
fontsize(22,"points");
grid on;

hold off;

f3 = figure(3);
f3_position = [p2X0+2*p2_scale*p2W,p2Y0,p2_scale*p2W,p2_scale*p2H];
f3.Position = f3_position;

%figure(3)
plot(t(2:Ntime)/sec_in_yr/1e6,angle_from_x(2:Ntime)); hold on;
xlabel("time (Myr)")
ylabel("angle between a and positive X-axis")
fontsize(22,"points");
grid on;

hold off;

f4 = figure(4);
f4_position = [p2X0,p2Y0+p2_scale*p2H,p2_scale*p2W,p2_scale*p2H];
f4.Position = f4_position;

%figure(4)
plot(t/sec_in_yr/1e6,finite_strain_q,'k--',linewidth=10); hold on;
plot(t/sec_in_yr/1e6,finite_strain_qps,'b+',markersize=30); hold on;
plot(t/sec_in_yr/1e6,finite_strain_qss,'rp',markersize=30); 
plot(t/sec_in_yr/1e6,finite_strain_qrr,'ko',markersize=30); 
legend("Lq",'Lqps','Lqss','Lqrr');
xlabel("time (Myr)")
ylabel("finite strain")
grid on;
fontsize(22,"points");
hold off;

f5 = figure(5);
f5_position = [p2X0+p2_scale*p2W,p2Y0+p2_scale*p2H,p2_scale*p2W,p2_scale*p2H];
f5.Position = f5_position;
%figure(5)
plot(t(2:Ntime)/sec_in_yr/1e6,angle_from_a(2:Ntime),'go',markersize=30); hold on;
plot(t(2:Ntime)/sec_in_yr/1e6,angle_from_a_qps(2:Ntime),'b+',markersize=30); 
plot(t(2:Ntime)/sec_in_yr/1e6,angle_from_a_qss(2:Ntime),'rp',markersize=30); 
plot(t(2:Ntime)/sec_in_yr/1e6,angle_from_a_qrr(2:Ntime),'ko',markersize=30); 
legend("Lq",'Lqps','Lqss','Lqrr');
xlabel("time (Myr)")
ylabel("angle between a and positive a-axis")
grid on;
fontsize(22,"points");
hold off;


f6 = figure(6);
f6_position = [p2X0+2*p2_scale*p2W,p2Y0+p2_scale*p2H,p2_scale*p2W,p2_scale*p2H];
f6.Position = f6_position;
plot(t/sec_in_yr/1e6,finite_strain_q-finite_strain_qps,'--',linewidth=3); hold on;
plot(t/sec_in_yr/1e6,finite_strain_qss,'r+',markersize=30); hold on;
legend("Lq-Lqps",'Lqss');
xlabel("time (Myr)")
ylabel("finite strain")
fontsize(22,"points");

grid on;
fprintf("sum of finite strain from Lqps: %.3f\n",sum(finite_strain_qps(1:Ntime-1)));
fprintf("delta finite strain end minus beginning: %.2f\n",(finite_strain(Ntime)-finite_strain(1)))
fprintf("percentage of pure shear contribution on finite strain: %.0d\n",sum(finite_strain_qps(1:Ntime-1))/(finite_strain(Ntime)-finite_strain(1))*100)
fprintf("rotation from rigid rotation in degrees%.3f\n",sum(angle_from_a_qrr(1:Ntime-1)))
fprintf("rotation from simple shear in degrees%.3f\n",sum(angle_from_a_qss(1:Ntime-1)))


%% 


    





