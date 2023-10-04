% Finite Strain Ellipse from different combination of strain history
%clear
clc
close all
clear
%% SIMPLE case
% S = .839;
% S = 0;
% W = -.1;
% W = 0;
% L = [0 S-W; S+W 0]; %McKenzie 1979 eq27 ; vel gradient tensor

epsilon_xy = 1e-15;% 1/s
dvxdy = epsilon_xy;
%L = [0 dvxdy; 0 0] % simple shear
%L = [0 0; dvxdy 0] % simple shear

%L = [0 dvxdy; -dvxdy 0] % rigid rotation counterclockwise 
L = [0 -dvxdy; dvxdy 0] % rigid rotation clockwise 

%L = [0 dvxdy; dvxdy 0] % L pure shear
%L = [dvxdy 0; 0 -dvxdy] % L compression vertically

E11 = L(1,1);
E22 = L(2,2);
E12 = 0.5*(L(1,2)+L(2,1));
E21 = E12;

R11 = 0;
R22 = 0;
R12 = 0.5*(L(2,1)-L(1,2));
R21 = -R12;

E = [E11 E12; E21 E22];
R = [R11 R12; R21 R22];

time_end = 30; % Myrs
sec_in_yr = 365.25*24*3600; %seconds
Ntime = 100;
t = linspace(1,time_end*1e6*sec_in_yr,Ntime);
dt = t(2)-t(1);
finite_strain = zeros(1,Ntime);
angle_from_x = zeros(1,Ntime);

Esum = E*1e6*sec_in_yr;
Rsum = R*1e6*sec_in_yr;

F = [1 0; 0 1]; % initial deformation gradient tensor
%F = [2 0; 0 1]; % initial deformation gradient tensor


U = (F*F'); %left stretch tensor
[V, D] = eig(U);
% need to sort result from eig because it could be random                                         
        % ref: https://www.mathworks.com/help/matlab/ref/eig.html                                         
[Dsortedvec,Dind] = sort(diag(D)); %ascending order                                               
Dsorted = D(Dind,Dind);
Vsorted = V(:,Dind);
% theta= atan(V(2,1)/V(1,1)); % in rad                                                            
theta= atan(Vsorted(2,2)/V(1,2)); % in rad; angle of the long axis   
FSEscale = 1; % make them look bigger
ellipse(FSEscale*sqrt(Dsorted(2,2)),FSEscale*sqrt(Dsortedvec(1,1)),theta,0, ...
0,'blue',100); 
axis equal;
%% set up the Lq for tracking sense of deformation relative to FSE
Lq = zeros(2,2,Ntime);
Lq_sum = zeros(2,2); %accummulated strain relative to a-b axis of FSE
L_sum = zeros(2,2); %accummulated strain relative to X-Y axis 

for i=1:Ntime
    i;
    F = F+dt*L*F
    U = (F*F'); %left stretch tensor
    [V, D] = eig(U);
    % need to sort result from eig because it could be random                                         
            % ref: https://www.mathworks.com/help/matlab/ref/eig.html                                         
    [Dsortedvec,Dind] = sort(diag(D)); %ascending order                                               
    Dsorted = D(Dind,Dind);
    Vsorted = V(:,Dind);
    theta= atan(Vsorted(2,2)/Vsorted(1,2)); % in rad; angle of the long axis   
    %% rotation matrix R [no, actually should use transformation matrix]
    % ref:https://www.continuummechanics.org/rotationmatrix.html#:~:text=Summary,is%20along%20for%20the%20ride.
    % ref: https://www.continuummechanics.org/coordxforms.html
    % ref: https://www.continuummechanics.org/transformmatrix.html
    %R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    Q = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    Lq(:,:,i) = Q*L*Q'; % rotate coordinate system from X-Y to a-b for strain rate/vel grad tensor to get Lq in a-b sys
    Lq_sum = Lq_sum + Lq(:,:,i)*dt;
    L_sum = L_sum + L*dt;

    theta/pi*180;
    FSEscale = 1; % make them look bigger
    disp("-------------------------------")
    disp("a*b=")
    sqrt(Dsortedvec(2))*sqrt(Dsortedvec(1))
    Vsorted(1,1)^2+Vsorted(2,1)^2
    disp("-------------------------------")

    ellipse(FSEscale*sqrt(Dsortedvec(2)),FSEscale*sqrt(Dsortedvec(1)),theta,0, ...
    0,'red',100); hold on;
    quiver(0,0,sqrt(Dsorted(2,2))*Vsorted(1,2),sqrt(Dsorted(2,2))*Vsorted(2,2),'red'); 
    quiver(0,0,sqrt(Dsorted(1,1))*Vsorted(1,1),sqrt(Dsorted(1,1))*Vsorted(2,1),'blue'); 

    axis equal;
    grid on;
    set(gca, 'YDir','reverse');

    finite_strain(i) = log(sqrt(Dsortedvec(2))/sqrt(Dsortedvec(1)));
    angle_from_x(i) = theta/pi*180;
    title(sprintf("Pure and Simple Finite strain=%.5f;long axis angle=%.5f;\ntime=%.1f[kyrs]",finite_strain(i),angle_from_x(i)),t(i)/sec_in_yr/1e3);
    pause(.05);
end

figure
plot(t/sec_in_yr/1e6,finite_strain); hold on;
xlabel("time (Myr)")
ylabel("finite strain log(a/b)")
figure
plot(t/sec_in_yr/1e6,angle_from_x); hold on;
xlabel("time (Myr)")
ylabel("angle between a and positive X-axis")

%% plotting the Lq_sum effect on original round ellipse
figure
F = [1 0; 0 1]; % initial deformation gradient tensor

U = (F*F'); %left stretch tensor
[V, D] = eig(U);
% need to sort result from eig because it could be random                                         
        % ref: https://www.mathworks.com/help/matlab/ref/eig.html                                         
[Dsortedvec,Dind] = sort(diag(D)); %ascending order                                               
Dsorted = D(Dind,Dind);
Vsorted = V(:,Dind);
% theta= atan(V(2,1)/V(1,1)); % in rad                                                            
theta= atan(Vsorted(2,2)/V(1,2)); % in rad; angle of the long axis   
FSEscale = 1; % make them look bigger
ellipse(FSEscale*sqrt(Dsorted(2,2)),FSEscale*sqrt(Dsortedvec(1,1)),theta,0, ...
0,'blue',100); 
axis equal;


F = [1 0; 0 1]; % initial deformation gradient tensor
%Wrong way to calculate the final FSE from below
% F = F + Lq_sum*F; %final ellipse from Lq_sum
% dF/dt = LF ---> Ft = exp(sum(L)dt) * F0
% the above is still wrong, you would need to do matrix differential eqs
% ref: https://en.wikipedia.org/wiki/Matrix_differential_equation
%% pure shear but compression in Y axis and extension in X axis
%F = [exp(1e-15*time_end*1e6*sec_in_yr) 0; 0 exp(-1e-15*time_end*1e6*sec_in_yr)] % for the L = [dvxdy 0; 0 -dvxdy]
F = [exp(1e-15*(time_end*1e6*sec_in_yr+dt)) 0; 0 exp(-1e-15*(time_end*1e6*sec_in_yr+dt))]; % for the L = [dvxdy 0; 0 -dvxdy]
%% pure shear 45 degree from XY
B = 1/1.4142;
A = -B;
Ap = 1/1.4142;
Bp = Ap;
Fps11 = A*exp(-dvxdy*(time_end*1e6*sec_in_yr+dt))*-0.7071 + B*exp(dvxdy*(time_end*1e6*sec_in_yr+dt))*0.7071;
Fps21 = A*exp(-dvxdy*(time_end*1e6*sec_in_yr+dt))* 0.7071 + B*exp(dvxdy*(time_end*1e6*sec_in_yr+dt))*0.7071;
Fps12 = Ap*exp(-dvxdy*(time_end*1e6*sec_in_yr+dt))*-0.7071 + Bp*exp(dvxdy*(time_end*1e6*sec_in_yr+dt))*0.7071;
Fps22 = Ap*exp(-dvxdy*(time_end*1e6*sec_in_yr+dt))* 0.7071 + Bp*exp(dvxdy*(time_end*1e6*sec_in_yr+dt))*0.7071;
%F = [Fps11 Fps12;Fps21 Fps22];

%% simple shear
Fsimpleshear21 = (time_end*1e6*sec_in_yr+dt)*dvxdy;
F = [1 0; Fsimpleshear21 1];

U = (F*F'); %left stretch tensor
[V, D] = eig(U);
% need to sort result from eig because it could be random                                         
        % ref: https://www.mathworks.com/help/matlab/ref/eig.html                                         
[Dsortedvec,Dind] = sort(diag(D)); %ascending order                                               
Dsorted = D(Dind,Dind);
Vsorted = V(:,Dind);
% theta= atan(V(2,1)/V(1,1)); % in rad                                                            
theta= atan(Vsorted(2,2)/V(1,2)); % in rad; angle of the long axis   
FSEscale = 1; % make them look bigger
ellipse(FSEscale*sqrt(Dsortedvec(2)),FSEscale*sqrt(Dsortedvec(1)),theta,0, ...
0,'red',100); 
axis equal;
grid on;
set(gca, 'YDir','reverse');


finite_strain_rela_FSE = log(sqrt(Dsortedvec(2))/sqrt(Dsortedvec(1)))
angle_rela_FSE_a = theta/pi*180


E11_FSE = Lq_sum(1,1);
E22_FSE = Lq_sum(2,2);
E12_FSE = 0.5*(Lq_sum(1,2)+Lq_sum(2,1));
E21_FSE = E12_FSE;

R11_FSE = 0;
R22_FSE = 0;
R12_FSE = 0.5*(Lq_sum(2,1)-Lq_sum(1,2));
R21_FSE = -R12_FSE;

E_FSE = [E11_FSE E12_FSE; E21_FSE E22_FSE]
R_FSE = [R11_FSE R12_FSE; R21_FSE R22_FSE]






