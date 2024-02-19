% SiStER_MAIN.m
%
% Simple Stokes solver with Exotic Rheologies
%
% Main routine doing initialization, time loop and outputs
%
%
% J.-A. Olive, B.Z. Klein, E. Mittelstaedt, M. Behn, G. Ito, S. Howell
% jaolive <at> ldeo.columbia.edu
% March 2011 - April 2017

close all
%<>
%Path_results=['/Users/tian_bc/repos/github/SiStER/dir_two_Layers/exp4/'];
%Path_results=['/Users/tian_bc/repos/github/SiStER/dir_profiling/'];
%Path_results=['/Users/tian_bc/repos/github/SiStER/dir_results_out/folder_HT_mylonite/input_extension1/']
Path_results=['/Users/tian_bc/repos/github/SiStER/dir_results_out/dir_Delamination_FSE/exp2_PI/'];
if exist(Path_results, 'dir')==0
    sprintf("path not existing")
    mkdir(Path_results)
end
addpath(Path_results)

% INITIALIZATION

% Input File: loads parameter values, model geometry, boundary conditions
% if exist('running_from_SiStER_RUN','var')==0
%      clear 
%      InpFil = input('Input file ? ','s');
%  end
%InpFil='SiStER_Input_File_extension_csdms.m'
InpFil='SiStER_Input_File_delamination_twoL.m'
run(InpFil)

% construct grid and initialize marker / node arrays
SiStER_Initialize

% BEGIN TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time=0;

for t=1:Nt % time loop
    
    disp(['STARTING ITERATION: ' num2str(t) ' out of ' num2str(Nt)])
    
    % update time
    time=time+dt_m;
    
    % Here we prepare nodal arrays to feed the Stokes solver 
    SiStER_material_props_on_nodes

    %%% SOLVE STOKES WITH NON-LINEAR RHEOLOGY HERE 
    SiStER_flow_solve
    
    % GET STRAIN RATE FROM CURRENT SOLUTION
    epsIIm=SiStER_interp_shear_nodes_to_markers(epsII_s,x,y,xm,ym,icn,jcn);
    
    % USE STRAIN RATE TO UPDATE STRESSES ON MARKERS
    SiStER_update_marker_stresses;

    % BUILD UP PLASTIC STRAIN IN YIELDING AREAS IF PLASTICITY IS ACTIVATED
    if (PARAMS.YNPlas==1) 
        SiStER_update_ep;
    end

    % update finite strain
    %before updating finite strain, saving the previous time step results for calc theta_PI
    ep1_xxm_old = ep1_xxm;
    ep1_xym_old = ep1_xym;
    ep1_yxm_old = ep1_yxm;
    ep1_yym_old = ep1_yym;
    
    SiStER_update_finite_strain;

    vxm_old = vxm;
    vym_old = vym;

    % MARKER ADVECTION, REMOVAL, AND ADDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    SiStER_move_remove_and_reseed_markers;
    % advect markers in current flow field
    % remove markers if necessary
    % add markers if necessary    
    SiStER_update_topography_markers;
    % here we do the same for the marker chain that keeps track of topography
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [PIm,theta_PIm,DthetaDt_m]=fabric_get_PI(epsIIm, ...
    vxm_old, vym_old, vxm,vym, ...
    ep1_xxm_old, ep1_xym_old, ep1_yxm_old, ep1_yym_old, ...
    ep1_xxm, ep1_xym, ep1_yxm, ep1_yym, dt_m);
  
    % OUTPUT VARIABLES OF INTEREST (prior to rotation & advection)
    if (mod(t,dt_out)==0 && dt_out>0) || t==1 || t==Nt % SAVING SELECTED OUTPUT
        disp('SAVING SELECTED VARIABLES TO OUTPUT FILE') 
        filename=num2str(t);
        [etam]=SiStER_interp_shear_nodes_to_markers(etas,x,y,xm,ym,icn,jcn); % to visualize viscosity on markers
        % saving finite strain tensor (pst)
        pst.ep1_xxm = ep1_xxm;
        pst.ep1_xym = ep1_xym;
        pst.ep1_yym = ep1_yym;
        pst.ep1_yxm = ep1_yxm;

        %pst.ep1_PIm = ep1_PIm;
        %save(filename,'X','Y','vx','vy','p','time','xm','ym','etam','rhom','rho','BC','etan','Tm','im','idm','epsIIm','sxxm','sxym','ep','epNH','icn','jcn','qd','topo_x','topo_y')

        %save([Path_results,filename],'X','Y','vx','vy','p','time','xm','ym','etam','rhom','rho','BC','etan','Tm','im','idm','epsIIm','sxxm','sxym','ep','epNH','icn','jcn','qd','topo_x','topo_y','pst')
        save([Path_results,filename],'X','Y','vx','vy','p','time','xm','ym','etam','rhom','rho','BC','etan','Tm','im','idm','epsIIm','sxxm','sxym','ep','epNH','icn','jcn','qd','topo_x','topo_y','pst','vxm','vym','dt_m','PIm','theta_PIm','DthetaDt_m')

    end
    
    % SET ADVECTION TIME STEP BASED ON CURRENT FLOW SOLUTION
    [dt_m]=SiStER_set_timestep(dx,dy,vx,vy,PARAMS);

    % ROTATE ELASTIC STRESSES IN CURRENT FLOW FIELD
    if (PARAMS.YNElast==1) 
        SiStER_rotate_stresses;
    end
    
    % EVOLVE TEMPERATURE FIELD THROUGH DIFFUSION
    if PARAMS.Tsolve==1
        SiStER_thermal_update;
    end

    %% MARKER ADVECTION, REMOVAL, AND ADDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %SiStER_move_remove_and_reseed_markers;
    %% advect markers in current flow field
    %% remove markers if necessary
    %% add markers if necessary    
    %SiStER_update_topography_markers
    %% here we do the same for the marker chain that keeps track of topography
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('---------------')
    disp(['END OF ITERATION: ' num2str(t) ' out of ' num2str(Nt) ' - SIMULATION TIME: ' num2str(time/365.25/24/3600/1000) ' kyrs.'])
    disp('--------------------------------')
    disp('--------------------------------')
    

end

disp('FIN')

    
