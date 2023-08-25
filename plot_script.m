close all;
% 
% %fastscatter(xm(im>1)/1e3,ym(im>1) /1e3,log10(epsIIm(im>1)),'markersize',2);

list_step = [1,10,20,30,40,60,80,100,120,140,160,180,200];
%list_step = [1];


sec_in_yr = 365*24*3600;

f1 = figure('Position',[0 0 600 400]);
f2 = figure('Position',[610 0 600 400]);
f3 = figure('Position',[0 460 600 400]);
f4 = figure('Position',[610 460 600 400]);
f5 = figure('Position',[1230 0 450 200]);
for i=1:length(list_step)
    i
    matFilename  = sprintf('%d.mat', list_step(i));
    fprintf(matFilename)
    matData      = load(fullfile(matFilename));
    %pcolor(matData.p)
    %pcolor(matData.rho)

    fprintf('\n time:%.1f kyrs \n',matData.time/sec_in_yr/1e3)

    %f1 = figure('Position',[0 0 600 400]);
    figure(f1)
    fastscatter(matData.xm/1e3,matData.ym/1e3,matData.rhom)
    axis ij;
    %axis equal;
    grid off;
    ylim([0 600]) 
    xlim([0 500])
    text_title = sprintf('Density[kg/m3] \n time:%.0f kyrs \n',matData.time/sec_in_yr/1e3)
    title(text_title);
    colorbar
    %movegui(f1,[0 0]);
    
    %f2 = figure('Position',[630 0 600 400]);
    figure(f2)
    fastscatter(matData.xm/1e3,matData.ym/1e3,matData.im,'markersize',2);
    set(gcf,'color','white')
    set(gca,'YDir','reverse')
    %axis equal
    colorbar
    text_title = sprintf('MATERIALs \n time:%.0f kyrs \n',matData.time/sec_in_yr/1e3)
    title(text_title);
    %title('MATERIALS \n')
    xlabel('X distance (km)','fontsize',15)
    ylabel('depth (km)','fontsize',15)
    ylim([0 600]) 
    xlim([0 500])
    %movegui(f2,[630 0]);

    %f3 = figure('Position',[0 420 600 400]);
    figure(f5)
    plot(matData.topo_x/1e3,matData.topo_y/1e3)
    xlabel('X distance (km)','fontsize',15)
    ylabel('topography (km)','fontsize',15)
    set(gca,'YDir','reverse')
    ylim([48 52]) 
    xlim([0 500])


    %f4 = figure('Position',[630 420 600 400]);
    figure(f4)
    fastscatter(matData.xm/1e3,matData.ym/1e3,matData.sxym,'markersize',2);
    xlabel('X distance (km)','fontsize',15)
    ylabel('depth (km)','fontsize',15)
    title('sxy')
    set(gca,'YDir','reverse')
    ylim([0 600]) 
    xlim([0 500])
    colorbar

    figure(f3)
    fastscatter(matData.xm/1e3,matData.ym/1e3,log10(matData.etam),'markersize',2);
    xlabel('X distance (km)','fontsize',15)
    ylabel('depth (km)','fontsize',15)
    title('etam')
    set(gca,'YDir','reverse')
    ylim([0 600]) 
    xlim([0 500])
    colorbar

    pause(0.3)
end

