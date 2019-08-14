function [ob]= ob_ini_fromTrasect(gr, trFname)
%[ob]= ob_ini_fromTrasect(gr, trFname)
% initialize observation structure 'ob' 
% from grid structure 'gr' and trasect file 'trFname' 
%INPUT/OUTPUT
% ob - observation datastructure 
% trFanme - finame for trasect, should be a build point file
% gr - grid strucutre with hgrid and vgrid fileds
%
% Sergey Frolov, June 2005

%initialise ob strucutre
ob 	= ob_ini(gr.hgrid, gr.vgrid);

if exist(trFname)~=2
    error(['file ' trFname ' doesnt exist'])
end

%read transect.bp file
[tr_bp]     = gr_readHGrid(trFname);

%assign xy points
ob.xy.x     = tr_bp.x;
ob.xy.y     = tr_bp.y;
ob.xy.flag  = 1;    %points loaded
%disp('transect loaded')

%compute interpolation weights
[ob]        = ob_comp_wxy(ob);
ob.xy.flag  = 2;    %weights computed
%disp('weights computed')

%assemble sparce observation matrix H
[ob]=ob_obxy2h(ob);
ob.xy.flag  = 3;    %matrix H assembled
%disp('matrix H assembled')

