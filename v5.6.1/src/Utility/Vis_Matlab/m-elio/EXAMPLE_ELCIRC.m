clear all

%%%%%%%%%%%%%%%%%%
%ELCIRC
%Data Format 3.0
%%%%%%%%%%%%%%%%%%

%links to files on CCALMR network, substitute appropriately
fn.runDir = '/home/workspace/ccalmr/hindcasts/2002-11-11/run/';

%%%%%%%%%%%%%%%%%%
% grid
%%%%%%%%%%%%%%%%%%
%read horiz. and vert. grids
gr.hgrid=gr_readHGrid(fullfile(fn.runDir,'hgrid.gr3'));
gr.vgrid=gr_readVGrid(fullfile(fn.runDir,'vgrid.in'));

%%%%%%%%%%%%%%%%%%
% binary files and graphics
%%%%%%%%%%%%%%%%%%

%reading, mapping and plotting 3d output
%velocity
 %read header of the elcirc file, need to do it once for each file type
h=eb_readHeader(fullfile(fn.runDir,'1_hvel.64'));
[d ts]=eb_readTimeStep(h,1); 
 %for 3D data map elcirc binary output to format [nPoint nLevs]
[u] = map_eb2hts(h,d(:,1));% map component u of velocity
gr_plot(gr.hgrid,u(:,end)); %surface velocity

%for elevation
h=eb_readHeader(fullfile(fn.runDir,'1_elev.61'));
[e ts]=eb_readTimeStep(h,1);
gr_plot(gr.hgrid,e); %elevation

%for a file different then h.fname
h.fname = fullfile(fn.runDir,'2_elev.61')
[e ts]=eb_readTimeStep(h,1);
gr_plot(gr.hgrid,e); %elevation


%%%%%%%%%%%%%%%%%%
%interpolation
%%%%%%%%%%%%%%%%%%

%triangulate the grid, since xy interpolation uses native matlab tools
%hence be careful with interpolation outside of the domain
gr=gr_tri(gr);
xy = [343400.0  288056.0];
[ob]= ob_ini_staFromPars(gr, xy);
% ob.xy.H is a sparse observation matrix

%horizontal interpolation gets a column of u values 
obs=ob.xy.H*u;
%vertical interpolation,
obs=interp1(gr.vgrid.zLevel-gr.vgrid.zMsl,obs,-2);
%can also use horizontal interpolation for vis. of transects


