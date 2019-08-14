function [ob]= ob_loadTrasect(ob, trFname)
%[ob]= ob_loadTrasect(ob, trFname)
% loads trasect file (*.bp) into ob datastructure
%INPUT/OUTPUT
% ob - observation datastructure as initialized by ob_ini
% trFanme - finame for trasect
% assumes that no data was loaded to "ob" yet
%
% Sergey Frolov, Dec 2004


if exist(trFname)~=2
    error(['file ' trFname ' doesnt exist'])
end

%read transect.bp file
[tr_bp]     = gr_readHGrid(trFname);

%assign xy points
ob.xy.x     = tr_bp.x;
ob.xy.y     = tr_bp.y;
ob.xy.flag  = 1;    %points loaded
disp('transect loaded')

%compute interpolation weights
[ob]        = ob_comp_wxy(ob);
ob.xy.flag  = 2;    %weights computed
disp('weights computed')

%assemble sparce observation matrix H
[ob]=ob_obxy2h(ob);
ob.xy.flag  = 3;    %matrix H assembled
disp('matrix H assembled')
