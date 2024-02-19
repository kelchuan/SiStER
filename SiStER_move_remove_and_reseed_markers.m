% SiStER Update Markers

%[xm_new,ym_new] = SiStER_advect_markers(x,y,xm,ym,dx,dy,dt_m,vx,vy);
[xm_new,ym_new, vxm_new,vym_new] = SiStER_advect_markers(x,y,xm,ym,dx,dy,dt_m,vx,vy); %  (Tian240219)  
xm=xm_new;
ym=ym_new;
vxm = vxm_new;  %(Tian240219)  
vym = vym_new;  % (Tian240219)  
% eliminate markers that left domain
Iin=find(xm<=xsize & xm>=0 & ym>=0 & ym<=ysize);

%msg2='  markers removed: ';
%msg=[msg2 num2str(length(xm)-length(Iin))];
%disp(msg)
xm=xm(Iin);
ym=ym(Iin);
im=im(Iin);
ep=ep(Iin);
epNH=epNH(Iin);
Tm=Tm(Iin);
idm=idm(Iin);
sxxm=sxxm(Iin);
sxym=sxym(Iin);
epsIIm=epsIIm(Iin);

% finite strain tensor                                                                                                                                                    
ep1_xxm = ep1_xxm(Iin) ;
ep1_yym = ep1_yym(Iin) ;
ep1_xym = ep1_xym(Iin) ;
ep1_yxm = ep1_yxm(Iin) ;

%velocity
vxm = vxm(Iin);
vym = vym(Iin);

% locate advected markers with respect to the eulerian grid
[quad,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);
    
% check for holes in the marker distribution, 
% patch with new markers if necessary
% those new markers immediately get assigned a value of phase (im), index 
% (idm) and accumulated plastic strain (ep), i.e., the 2 variables that never get
% passed to nodes. 
%[xm, ym, im, Ifix, mp, ep, idm, Tm, sxxm, sxym, epNH, epsIIm, ep1_xxm,ep1_yym,ep1_xym,ep1_yxm]=SiStER_patch_marker_holes(icn,jcn,quad,Nx,Ny,Mquad,Mquad_crit,xm,ym,x,y,dx,dy,im,ep,idm,Tm,sxxm,sxym,epNH, epsIIm, ep1_xxm,ep1_yym,ep1_xym,ep1_yxm);
[xm, ym, im, Ifix, mp, ep, idm, Tm, sxxm, sxym, epNH, epsIIm, ep1_xxm,ep1_yym,ep1_xym,ep1_yxm,vxm,vym]=SiStER_patch_marker_holes(icn,jcn,quad,Nx,Ny,Mquad,Mquad_crit,xm,ym,x,y,dx,dy,im,ep,idm,Tm,sxxm,sxym,epNH, epsIIm, ep1_xxm,ep1_yym,ep1_xym,ep1_yxm,vxm,vym);    

% then they get assigned P, epsII and stresses from grid values

if min(Ifix)>0
    

    xmFIX=xm(Ifix);
    ymFIX=ym(Ifix);
    
    % pass temperature, pressure, strain rate and stresses to the new
    % markers from their nodal values
    % locate new markers with respect to the eulerian grid
    [quadFIX,icnFIX,jcnFIX] = SiStER_locate_markers_in_grid(xmFIX,ymFIX,x,y,dx,dy);
    
    
    [temp]=SiStER_interp_normal_nodes_to_markers(p,xc,yc,xmFIX,ymFIX,icnFIX,jcnFIX);
    pm(Ifix)=temp; % pressure
    
    
end
    
    
% locate all markers with respect to the eulerian grid
[qd,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);


