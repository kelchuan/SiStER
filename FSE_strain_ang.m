function [strain,angle] = FSE_strain_ang(matData, markeri, scale, color)
    %markeri
    F11 = matData.pst.ep1_xxm(markeri);
    F12 = matData.pst.ep1_xym(markeri);
    F21 = matData.pst.ep1_yxm(markeri);
    F22 = matData.pst.ep1_yym(markeri);
    FiniteStrain = [F11 F12; F21 F22];
    %left_stretch_tensor =  (FiniteStrain*FiniteStrain')^0.5;
    % left stretch tensor Uik = (Fij*Fjk')^(1/2)
    % take the s
    left_stretch_tensor =  (FiniteStrain*FiniteStrain');%sqrt for only D
    [V, D] = eig(left_stretch_tensor);
    % need to sort result from eig because it could be random
    % ref: https://www.mathworks.com/help/matlab/ref/eig.html
    [Dsortedvec,Dind] = sort(diag(D)); %ascending order
    Dsorted = D(Dind,Dind);
    Vsorted = V(:,Dind);
    % theta= atan(V(2,1)/V(1,1)); % in rad
    %theta= atan(Vsorted(2,2)/V(1,2)); % in rad; angle of the long axis
    theta= atan(Vsorted(2,2)/Vsorted(1,2)); % in rad; angle of the long axis

    FSEscale = scale; % make them look bigger
    ellipse(FSEscale*sqrt(Dsortedvec(2)),FSEscale*sqrt(Dsortedvec(1)),theta,matData.xm(markeri),matData.ym(markeri),color);
    
    [u,v] = pol2cart(theta,FSEscale*sqrt(Dsortedvec(2)));
    quiver(matData.xm(markeri),matData.ym(markeri),u,v,'red')
    hold on;
    strain = log(sqrt(Dsortedvec(2))/sqrt(Dsortedvec(1)));
    angle = theta/pi*180;
    %title(sprintf("Pure and Simple Finite strain=%.5f; long axis angle=%.5f",strain,angle));
end