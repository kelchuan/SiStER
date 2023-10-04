% Finite Strain Ellipse from different combination of strain history
%clear
%clc
%close all

%% SIMPLE case
S = .839;
S = 0;
W = -.1;
W = 0;
L = [0 S-W; S+W 0]; %McKenzie 1979 eq27 ; vel gradient tensor
L = [0 0;.839 0]



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




F = F+L*F;

U = (F*F'); %left stretch tensor
[V, D] = eig(U);
% need to sort result from eig because it could be random                                         
        % ref: https://www.mathworks.com/help/matlab/ref/eig.html                                         
[Dsortedvec,Dind] = sort(diag(D)); %ascending order                                               
Dsorted = D(Dind,Dind);
Vsorted = V(:,Dind);
% theta= atan(V(2,1)/V(1,1)); % in rad                                                            
theta= atan(Vsorted(2,2)/V(1,2)); % in rad; angle of the long axis   
%FSEscale = 1; % make them look bigger
ellipse(FSEscale*sqrt(Dsorted(2,2)),FSEscale*sqrt(Dsortedvec(1,1)),theta,0, ...
0,'red',100);
axis equal;
grid on;

finite_strain = log(sqrt(Dsorted(2,2))/sqrt(Dsortedvec(1,1)));
angle_from_x = theta/pi*180;
title(sprintf("Pure and Simple Finite strain=%.5f; long axis angle=%.5f",finite_strain,angle_from_x));

