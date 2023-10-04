%plot Finite Strain Ellipse from SiStER
close all;
%list_step = [1,10,20,30,40,60,80,100,120,140,160,180,200];
max_step = 200;

Xmin = 0;
Xmax = 500;
Ymin = 0;
Ymax = 600;

sec_in_yr = 365*24*3600;
f1 = figure('Position',[0 0 1000 1200]);
% f1 = figure('Position',[0 0 600 400]);
% f2 = figure('Position',[610 0 600 400]);
% f3 = figure('Position',[0 460 600 400]);
% f4 = figure('Position',[610 460 600 400]);
% f5 = figure('Position',[1230 0 450 200]);
%for i=1:length(list_step)


for i=10:10:max_step
    i
    %matFilename  = sprintf('%d.mat', list_step(i));
    matFilename  = sprintf('%d.mat', i);
    fprintf(matFilename)
    matData      = load(fullfile(matFilename));
    %pcolor(matData.p)
    %pcolor(matData.rho)
    fprintf('\n time:%.1f kyrs \n',matData.time/sec_in_yr/1e3)
    
    
    % to find the index of marker of interest
    % first use brush to select the marker and output their XY
    % then find(abs(xm-74.9523*1e3) < .1) and check ym(ind) is consistent
    % then use that ind below for calculate FSE

    figure(f1)
    %fastscatter(matData.xm(im>1)/1e3,matData.ym(im>1)/1e3,log10(matData.epsIIm(im>1)),'markersize',2);
    %fastscatter(matData.xm/1e3,matData.ym/1e3,log10(matData.epsIIm),'markersize',2);
    %fastscatter(matData.xm/1e3,matData.ym/1e3,matData.rhom,'markersize',2);
    fastscatter(matData.xm/1e3,matData.ym/1e3,matData.im,'markersize',12);
    %pcolor(matData.X/1e3, matData.Y/1e3, matData.rho)
    set(gcf,'color','white')
    set(gca,'YDir','reverse')
    axis equal
    colorbar
    %text_title = sprintf('SECOND INVARIANT OF STRAIN RATE (s^{-1}) \n time:%.0f kyrs \n',matData.time/sec_in_yr/1e3);
    text_title = sprintf('MATERIALs \n time:%.0f kyrs \n',matData.time/sec_in_yr/1e3);
    title(text_title);
    xlabel('cross-axis distance (km)','fontsize',15)
    ylabel('depth (km)','fontsize',15)

    %prep for FSE
    for markeri=1:20:length(matData.xm)
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
        theta= atan(Vsorted(2,2)/V(1,2)); % in rad; angle of the long axis
        FSEscale = 6; % make them look bigger
        %ellipse(FSEscale*sqrt(D(1,1)),FSEscale*sqrt(D(2,2)),theta,matData.xm(markeri)/1e3,matData.ym(markeri)/1e3,'white');
        ellipse(FSEscale*sqrt(Dsorted(2,2)),FSEscale*sqrt(Dsortedvec(1,1)),theta,matData.xm(markeri)/1e3,matData.ym(markeri)/1e3,'white');
        hold on;
    end

    ylim([Ymin Ymax]) 
    xlim([Xmin Xmax])

    pause(0.3)
end
