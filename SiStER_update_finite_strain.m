 
%=========================================================================
% PLASTIC (EXCLUDING ELASTIC) STRAIN ACCUMULATION
%% ---> FINITE STRAIN update
% Not sure whether to use the actual non-elastic strain (i.e., with
% current stresses) or the approximate, in which stress is assumed to equal
% yield stress
% G.Ito 8/16; JAO 9/15 for non-healed ep, fixed by JAO 4/2017
% X.Tian 8/23-- modified for finite strain 
%=========================================================================


%[psim]=SiStER_get_dilatancy(im,ep,MAT);
%dep_s=zeros(size(epsII_s));
%dep_s(s_nodes_yield) = dt_m.*max(epsII_s(s_nodes_yield)-...
%(yield_s(s_nodes_yield)-sqrt(sxxOLD_s(s_nodes_yield).^2+sxyOLD(s_nodes_yield).^2))./(2.*Gs(s_nodes_yield).*dt_m),min(epsII_s(:))*1e-6);
%[psrm]=SiStER_interp_shear_nodes_to_markers(reshape(dep_s,Ny,Nx),x,y,xm,ym,icn,jcn);
%depm = psrm / dt_m; % plastic strain rate (second invariant of deviatoric plastic strain rate tensor)

[rrs]=SiStER_get_rotation_rate(vx,vy,dx,dy,BC); % full rotation rate on shear nodes
[om]=SiStER_interp_shear_nodes_to_markers(rrs,x,y,xm,ym,icn,jcn); % rotation rate on markers
[EXYm]=SiStER_interp_shear_nodes_to_markers(EXY,x,y,xm,ym,icn,jcn); % EXY rate on markers
[EXXm]=SiStER_interp_normal_nodes_to_markers(EXX,xc,yc,xm,ym,icn,jcn); % EXX rate on markers

%% VELOCITY GRADIENT TENSOR 
% L = [∂vx/∂x   ∂vx/dy; ∂vy/dx - ∂vx/∂x]
% ∂vx/∂x = EXX; ∂vy/∂y = - ∂vx/∂x % continuity
% EXY = .5 (∂vx/∂y + ∂vy/dx)
% rotation rate om = .5 (∂vx/∂y - ∂vy/dx)
% rotation rate om = - .5 (∂vx/∂y - ∂vy/dx)
%<bug> om(i,j)=(-dvxdy + dvydx)/2; (clockwise positive; see SiStER_get_rotation_rate.m )
%  ∂vx/dy = EXY+om;    ∂vy/dx = ( EXY-om);

dep1_xxm = EXXm;
%dep1_xym = EXYm+om; %<bug>
%dep1_yxm = EXYm-om; %<bug>
dep1_xym = EXYm-om;
dep1_yxm = EXYm+om;
dep1_yym = -EXXm; 

%% GET PRINCIPAL STRESS ORIENTATIONS
% sIIm = sqrt(sxxm.^2+sxym.^2);
% ST = (sIIm + sxxm)./(sxym);
% v1 = ones(size(ST));
% v2 = ones(size(ST));
% v1(sxym~=0) = 1./sqrt(1+ST(sxym~=0).^2);
% v2(sxym~=0) = -ST(sxym~=0)./sqrt(1+ST(sxym~=0).^2);
% v1(sxym==0 & sxxm >= 0) = 0;
% v2(sxym==0 & sxxm >= 0) = 1;
% v1(sxym==0 & sxxm < 0) = 1;
% v2(sxym==0 & sxxm < 0) = 0;

%% VELOCITY GRADIENT TENSOR INCREMENTS
% % mode-I plastic strain
% dep1_xxm = 2*depm.*sind(psim).*v2.^2;
% dep1_yym = 2*depm.*sind(psim).*v1.^2;
% dep1_xym = -2*depm.*sind(psim).*v2.*v1 + om;
% dep1_yxm = -2*depm.*sind(psim).*v2.*v1 - om;
% % mode-II plastic strain
% gofpsi = (1-(4/3)*sind(psim).^2)./(sind(psim)+sqrt(1-(1/3)*sind(psim).^2));
% dep2_xxm = gofpsi.*depm.*(v2.^2-v1.^2);
% dep2_yym = gofpsi.*depm.*(v1.^2-v2.^2);
% dep2_xym = -2*gofpsi.*depm.*v1.*v2 + om;
% dep2_yxm = -2*gofpsi.*depm.*v1.*v2 - om;



%% ACTUAL FINITE STRAIN UPDATE
ep1_xxm_new = ep1_xxm + dt_m*(dep1_xxm.*ep1_xxm + dep1_xym.*ep1_yxm);
ep1_xym_new = ep1_xym + dt_m*(dep1_xxm.*ep1_xym + dep1_xym.*ep1_yym);
ep1_yxm_new = ep1_yxm + dt_m*(dep1_yxm.*ep1_xxm + dep1_yym.*ep1_yxm);
ep1_yym_new = ep1_yym + dt_m*(dep1_yxm.*ep1_xym + dep1_yym.*ep1_yym);

ep1_xxm = ep1_xxm_new; 
ep1_xym = ep1_xym_new;
ep1_yxm = ep1_yxm_new;
ep1_yym = ep1_yym_new;

%ep2_xxm_new = ep2_xxm + dt_m*(dep2_xxm.*ep2_xxm + dep2_xym.*ep2_yxm);
%ep2_xym_new = ep2_xym + dt_m*(dep2_xxm.*ep2_xym + dep2_xym.*ep2_yym);
%ep2_yym_new = ep2_yym + dt_m*(dep2_yxm.*ep2_xym + dep2_yym.*ep2_yym);
%ep2_yxm_new = ep2_yxm + dt_m*(dep2_yxm.*ep2_xxm + dep2_yym.*ep2_yxm);
%ep2_xxm = ep2_xxm_new; 
%ep2_xym = ep2_xym_new;
%ep2_yym = ep2_yym_new;
%ep2_yxm = ep2_yxm_new;



disp('** finite strain tensor UPDATED **')





