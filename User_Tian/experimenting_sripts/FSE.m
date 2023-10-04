% FSE.m
% '''given deformation gradient tensor F, 
%  calculate long axis a,short axis b, 
%  eigen vectors eigV, 
%  angle between a and X-axis theta (in rad), 
%  transformation matrix Q based on theta
%  of a FSE '''
function [a,b,eigV,theta,Q]=FSE(F)  
    %U = (F*F'); %left stretch tensor
    U = (F'*F); %right stretch tensor
    [V, D] = eig(U);
    % need to sort result from eig because it could be random                                         
    % ref: https://www.mathworks.com/help/matlab/ref/eig.html                                         
    [Dsortedvec,Dind] = sort(diag(D)); %ascending order                                               
    Dsorted = D(Dind,Dind);
    Vsorted = V(:,Dind);
    theta= atan(Vsorted(2,2)/Vsorted(1,2)); % in rad; angle of the long axis relative to X-axis 
    % rotation matrix R [no, actually should use transformation matrix]
    % ref:https://www.continuummechanics.org/rotationmatrix.html#:~:text=Summary,is%20along%20for%20the%20ride.
    % ref: https://www.continuummechanics.org/coordxforms.html
    % ref: https://www.continuummechanics.org/transformmatrix.html
    %R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    Q = [cos(theta) sin(theta); -sin(theta) cos(theta)]; %transformation matrix
    a = sqrt(Dsortedvec(2)); % long axis of FSE
    b = sqrt(Dsortedvec(1)); % short axis of FSE
    eigV = Vsorted;          % 2 by 2 matrix of eigenvectors of a and b
end