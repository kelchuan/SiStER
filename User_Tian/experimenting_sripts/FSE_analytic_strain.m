% Finite Strain Ellipse from different combination of strain history
%clear
%clc
close all

%% SIMPLE case
% S = .839;
% S = 0;
% W = -.1;
% W = 0;
% L = [0 S-W; S+W 0]; %McKenzie 1979 eq27 ; vel gradient tensor

epsilon_xy = 1e-15;% 1/s
dvxdy = epsilon_xy;
L = [0 dvxdy; 0 0] % L per 
L = [0 0; dvxdy 0] % L per 

time_end = 30; % Myrs
sec_in_yr = 365.25*24*3600; %seconds
Ntime = 100;
t = linspace(1,time_end*1e6*sec_in_yr,Ntime);
dt = t(2)-t(1);
finite_strain = zeros(1,Ntime);
angle_from_x = zeros(1,Ntime);

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
0,'blue',100); hold on;
axis equal;

for i=1:Ntime
    i
    F = F+dt*L*F;
    U = (F*F'); %left stretch tensor
    [V, D] = eig(U);
    % need to sort result from eig because it could be random                                         
            % ref: https://www.mathworks.com/help/matlab/ref/eig.html                                         
    [Dsortedvec,Dind] = sort(diag(D)); %ascending order                                               
    Dsorted = D(Dind,Dind);
    Vsorted = V(:,Dind);
    % theta= atan(V(2,1)/V(1,1)); % in rad                                                            
    %theta= atan(Vsorted(2,2)/Vsorted(1,2)); % in rad; angle of the long axis   
    theta= atan(Vsorted(2,2)/Vsorted(1,2)); % in rad; angle of the long axis   

    theta/pi*180
    FSEscale = 1; % make them look bigger
    ellipse(FSEscale*sqrt(Dsorted(2,2)),FSEscale*sqrt(Dsortedvec(1,1)),theta,0, ...
    0,'red',100); hold on;
    quiver(0,0,sqrt(Dsorted(2,2))*Vsorted(1,2),sqrt(Dsorted(2,2))*Vsorted(2,2),'red'); 
    quiver(0,0,sqrt(Dsorted(1,1))*Vsorted(1,1),sqrt(Dsorted(1,1))*Vsorted(2,1),'blue'); 

    axis equal;
    grid on;
    set(gca, 'YDir','reverse');

    finite_strain(i) = log(sqrt(Dsorted(2,2))/sqrt(Dsortedvec(1,1)));
    angle_from_x(i) = theta/pi*180;
    title(sprintf("Pure and Simple Finite strain=%.5f; long axis angle=%.5f",finite_strain(i),angle_from_x(i)));
    pause(.05)
end

figure
plot(t/sec_in_yr/1e6,finite_strain); hold on;
xlabel("time (Myr)")
ylabel("finite strain log(a/b)")
figure
plot(t/sec_in_yr/1e6,angle_from_x); hold on;
xlabel("time (Myr)")
ylabel("angle between a and positive X-axis")


