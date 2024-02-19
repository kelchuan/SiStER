function [PIm,theta_PIm,DthetaDt_m]=fabric_get_PI(epsIIm, ...
    vxm_old, vym_old, vxm,vym, ...
    ep1_xxm_old, ep1_xym_old, ep1_yxm_old, ep1_yym_old, ...
    ep1_xxm, ep1_xym, ep1_yxm, ep1_yym, dt_m)
% Kaminski and Ribe 2002 G3 equation 9

%% for PI factor  %equation 9 of Kaminski&Ribe-2002-G3-Timescales ...
%PIm = 1/epsilon_dot * (Dtheta/Dt)
PIm = zeros(size(epsIIm));
theta_PIm = zeros(size(epsIIm)); % not sure if it is degrees or radian, for now assume degrees
theta_PIm_old = zeros(size(epsIIm));

%% [epsilon_dot is the largest eigenvalue of strain rate tensor --> epsIIm]
one_over_epsilon_dotm = (ones(size(epsIIm))./epsIIm);

%% NEW velocity and FSE
%% get unit vectors of velocity fields
velocity_magnitude = sqrt(vxm.^2 + vym.^2);
univec_vxm = vxm./velocity_magnitude;
univec_vym = vym./velocity_magnitude;

%% get unit vectors of long axis of FSE
FSE_laxis_vxm = zeros(size(epsIIm));
FSE_laxis_vym = zeros(size(epsIIm));

for markeri=1:length(epsIIm)
    %markeri
    F11 = ep1_xxm(markeri);
    F12 = ep1_xym(markeri);
    F21 = ep1_yxm(markeri);
    F22 = ep1_yym(markeri);
    FiniteStrain = [F11 F12; F21 F22];
    %left_stretch_tensor =  (FiniteStrain*FiniteStrain')^0.5;
    % left stretch tensor Uik = (Fij*Fjk')^(1/2)
    left_stretch_tensor =  (FiniteStrain*FiniteStrain');%sqrt for only D
    [V, D] = eig(left_stretch_tensor);
    % need to sort result from eig because it could be random
    % ref: https://www.mathworks.com/help/matlab/ref/eig.html
    [Dsortedvec,Dind] = sort(diag(D)); %ascending order
    Dsorted = D(Dind,Dind);
    Vsorted = V(:,Dind);
    % theta= atan(V(2,1)/V(1,1)); % in rad
    theta= atan(Vsorted(2,2)/Vsorted(1,2)); % in rad; angle of the long axis
    %Vsorted(2,2)^2 + Vsorted(1,2)^2
    FSE_laxis_vxm(markeri) = Vsorted(1,2);
    FSE_laxis_vym(markeri) = Vsorted(2,2);
end

%% OLD velocity and FSE
%% get unit vectors of velocity fields
velocity_magnitude_old = sqrt(vxm_old.^2 + vym_old.^2);
univec_vxm_old = vxm_old./velocity_magnitude_old;
univec_vym_old = vym_old./velocity_magnitude_old;

%% get unit vectors of long axis of FSE
FSE_laxis_vxm_old = zeros(size(ep1_xxm_old));
FSE_laxis_vym_old = zeros(size(ep1_xxm_old));

for markeri=1:length(ep1_xxm_old)
    %markeri
    F11 = ep1_xxm_old(markeri);
    F12 = ep1_xym_old(markeri);
    F21 = ep1_yxm_old(markeri);
    F22 = ep1_yym_old(markeri);
    FiniteStrain = [F11 F12; F21 F22];
    %left_stretch_tensor =  (FiniteStrain*FiniteStrain')^0.5;
    % left stretch tensor Uik = (Fij*Fjk')^(1/2)
    left_stretch_tensor =  (FiniteStrain*FiniteStrain');%sqrt for only D
    [V, D] = eig(left_stretch_tensor);
    % need to sort result from eig because it could be random
    % ref: https://www.mathworks.com/help/matlab/ref/eig.html
    [Dsortedvec,Dind] = sort(diag(D)); %ascending order
    Dsorted = D(Dind,Dind);
    Vsorted = V(:,Dind);
    % theta= atan(V(2,1)/V(1,1)); % in rad
    theta= atan(Vsorted(2,2)/Vsorted(1,2)); % in rad; angle of the long axis
    %Vsorted(2,2)^2 + Vsorted(1,2)^2
    FSE_laxis_vxm_old(markeri) = Vsorted(1,2);
    FSE_laxis_vym_old(markeri) = Vsorted(2,2);
end

%% calculate theta_PI
for markeri=1:length(ep1_xxm_old)
    vector_dot = univec_vxm(markeri)*FSE_laxis_vxm(markeri)...
        +univec_vym(markeri)*FSE_laxis_vym(markeri);
    theta_PIm(markeri) = acosd(vector_dot);

    vector_dot_old = univec_vxm_old(markeri)*FSE_laxis_vxm_old(markeri)...
        +univec_vym_old(markeri)*FSE_laxis_vym_old(markeri);
    theta_PIm_old(markeri) = acosd(vector_dot_old);
end
DthetaDt_m = abs((theta_PIm - theta_PIm_old) ./ dt_m);
PIm = one_over_epsilon_dotm.*DthetaDt_m;



