clear all

%%%%%%%%%%%%%%%%%%
%SELFE, Sigma-Z grid
%Data Format 5.0
%%%%%%%%%%%%%%%%%%

%links to files on CCALMR network, substitute appropriately
fn.runDir = '/home/workspace/ccalmr/hindcasts/2004-24-14/run/';

%%%%%%%%%%%%%%%%%%
% grid
%%%%%%%%%%%%%%%%%%
%read horiz. grid from a double precision hgrid.gr3
gr.hgrid=gr_readHGrid(fullfile(fn.runDir,'hgrid.gr3'));
%vertical grid
%gr_readVGrid doesn't work on sigma grids, instead use file header
h=sz_readHeader(fullfile(fn.runDir,'1_elev.61'));
gr.vgrid=h.vgrid;

%%%%%%%%%%%%%%%%%%
% binary files and graphics
%%%%%%%%%%%%%%%%%%

%reading, mapping and plotting 3d output
%velocity
 %read header of the elcirc file, need to do it once for each file type
h=sz_readHeader(fullfile(fn.runDir,'1_salt.63'));
[d ts]=sz_readTimeStep(h,1); 
 %for 3D data map elcirc binary output to format [nPoint nLevs]
[u] = map_sz2hts(h,d(:),1);% map component u of velocity
gr_plot(gr.hgrid,u(:,end)); %surface velocity

%for elevation
h=sz_readHeader(fullfile(fn.runDir,'1_elev.61'));
[e ts]=sz_readTimeStep(h,1);
gr_plot(gr.hgrid,e); %elevation

%for a file different then h.fname
[e ts]=sz_readTimeStep(h,1,fullfile(fn.runDir,'2_elev.61'));
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
obs=ob.xy.H*double(u);
%vertical interpolation,
%compute z levels first
[z] = sz_computeZlevels(ob.xy.H*gr.hgrid.depth, ob.xy.H*double(e), gr.vgrid);
obs=interp1(z,obs,-2)
%can also use horizontal interpolation for vis. of transects


