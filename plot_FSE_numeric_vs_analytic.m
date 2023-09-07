%plot Finite Strain Ellipse from SiStER
close all;
clc;
clear;
%list_step = [1,10,20,30,40,60,80,100,120,140,160,180,200];
max_step = 100;

Xmin = 0;
Xmax = 15;
Ymin = 0;
Ymax = 3;

sec_in_yr = 365*24*3600;
f1 = figure('Position',[0 0 1600 800]);

% list and sort all .mat data file
listmat = dir(fullfile('*.mat'));
listmatname = {listmat.name};
listmatname = natsortfiles(listmatname); %sort name so 100.mat is not before 20.mat
N_steps = length(listmatname); % N time step to plot

% for benchmarking, choose limited amount of ellipses
% for having a matData0 for init case
matData0     = load(fullfile(listmatname{1}));
N_ellipse = 50;
Ind_ellipse = linspace(1,length(matData0.xm),N_ellipse);
data_FSE(N_steps ,N_ellipse) = struct();


for i=1:N_steps
    i
    matFilename  = listmatname{i};
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
    
    %fastscatter(matData.xm,matData.ym,matData.im,'markersize',12);
    pcolor(matData.X, matData.Y, matData.vx); hold on;
    set(gcf,'color','white')
    set(gca,'YDir','reverse')
    axis equal
    colorbar
    %text_title = sprintf('SECOND INVARIANT OF STRAIN RATE (s^{-1}) \n time:%.0f kyrs \n',matData.time/sec_in_yr/1e3);
    text_title = sprintf('MATERIALs \n time:%.0f kyrs \n',matData.time/sec_in_yr/1e3);
    title(text_title);
    xlabel('cross-axis distance (m)','fontsize',15)
    ylabel('depth (m)','fontsize',15)

    %prep for FSE
    for markeri=1:2e3:length(matData0.xm)
    %for markeri=5e4
        % %markeri
        % F11 = matData.pst.ep1_xxm(markeri);
        % F12 = matData.pst.ep1_xym(markeri);
        % F21 = matData.pst.ep1_yxm(markeri);
        % F22 = matData.pst.ep1_yym(markeri);
        % FiniteStrain = [F11 F12; F21 F22];
        % %left_stretch_tensor =  (FiniteStrain*FiniteStrain')^0.5;
        % % left stretch tensor Uik = (Fij*Fjk')^(1/2)
        % % take the s
        % left_stretch_tensor =  (FiniteStrain*FiniteStrain');%sqrt for only D
        % [V, D] = eig(left_stretch_tensor);
        % % need to sort result from eig because it could be random
        % % ref: https://www.mathworks.com/help/matlab/ref/eig.html
        % [Dsortedvec,Dind] = sort(diag(D)); %ascending order
        % Dsorted = D(Dind,Dind);
        % Vsorted = V(:,Dind);
        % % theta= atan(V(2,1)/V(1,1)); % in rad
        % theta= atan(Vsorted(2,2)/V(1,2)); % in rad; angle of the long axis
        % FSEscale = .3; % make them look bigger
        % ellipse(FSEscale*sqrt(Dsorted(2,2)),FSEscale*sqrt(Dsortedvec(1,1)),theta,matData.xm(markeri),matData.ym(markeri),'red');
        % hold on;
        [strain,angle] = FSE_strain_ang(matData, markeri, 0.2, 'red')
    end
    ylim([Ymin Ymax]) 
    xlim([Xmin Xmax])
    pause(0.3)
    
    for j=1:N_ellipse
        i
        j
        data_FSE(i,j).time = matData.time;
        data_FSE(i,j).index = floor(Ind_ellipse(j));
        [FSEstrain,FSEangle]=FSE_strain_ang(matData, data_FSE(i,j).index, 0, 'red');
        data_FSE(i,j).strain = FSEstrain;
        data_FSE(i,j).angle = FSEangle;
    end
    numeric_ind_test = 10; %<>
    FSE_strain_ang(matData, data_FSE(i,numeric_ind_test).index, 1, 'blue');
end





%% analytic solution with the numerics
epsilon_xy = 1e-15;% 1/s
dvxdy = epsilon_xy;
L = [0 dvxdy; 0 0] % L per 
time_end = 30; % Myrs
sec_in_yr = 365.25*24*3600; %seconds
Ntime = 100;
t = linspace(1,time_end*1e6*sec_in_yr,Ntime);
dt = t(2)-t(1);
finite_strain = zeros(1,Ntime);
angle_from_x = zeros(1,Ntime);

F = [1 0; 0 1]; % initial deformation gradient tensor
for i=1:Ntime
    F = F+dt*L*F;
    U = (F*F'); %left stretch tensor
    [V, D] = eig(U);
    % need to sort result from eig because it could be random                                         
            % ref: https://www.mathworks.com/help/matlab/ref/eig.html                                         
    [Dsortedvec,Dind] = sort(diag(D)); %ascending order                                               
    Dsorted = D(Dind,Dind);
    Vsorted = V(:,Dind);
    theta= atan(Vsorted(2,2)/V(1,2)); % in rad; angle of the long axis   
    finite_strain(i) = log(sqrt(Dsorted(2,2))/sqrt(Dsortedvec(1,1)));
    angle_from_x(i) = theta/pi*180;
    %title(sprintf("Pure and Simple Finite strain=%.5f; long axis angle=%.5f",finite_strain(i),angle_from_x(i)));
end


%% numerical part
%data_FSE(1,numeric_ind_test).index
numeric_time = cell2mat({data_FSE(:,numeric_ind_test).time}');
numeric_strain = cell2mat({data_FSE(:,numeric_ind_test).strain}');
numeric_angle = cell2mat({data_FSE(:,numeric_ind_test).angle}')



figure
plot(t/sec_in_yr/1e6,finite_strain); hold on;
plot(numeric_time/sec_in_yr/1e6, numeric_strain,'b+',markersize=20); hold on;
xlabel("time (Myr)")
ylabel("finite strain log(a/b)")



figure
plot(t/sec_in_yr/1e6,angle_from_x); hold on;
plot(numeric_time/sec_in_yr/1e6, 90-numeric_angle,'b+',markersize=20); hold on;
%<?> not sure about why the angle needs 90- yet!
xlabel("time (Myr)")
ylabel("angle between a and positive X-axis")




