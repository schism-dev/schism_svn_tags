function [ ob ] = ob_ini( hgrid, vgrid )
%[ ob ] = ob_ini( hgrid, vgrid )
% initializes an ob datastructure 
% that defines the mapping (observation strategy) 
% from data defined on hgrid/vgrid to stations/transects ets defined on .xy .z and .xyz grids
%INPUT/OUTPUT
% hgrid as in gr.[hgrid|ecenters|sideCent] 
% vgrid as in gr.vgrid
% assumes that hgrid is triangulated
%
% Sergey Frolov, Dec 2004

ob.gr.hgrid = hgrid;
ob.gr.vgrid = vgrid;

try     %check that hgrid has triangulation available
    size(ob.gr.hgrid.tri);
catch   %if no tri is available the error is thrown, and we catch it and assign tri
    tri         = delaunay(ob.gr.hgrid.nodes(:,2),ob.gr.hgrid.nodes(:,3));
    ob.gr.hgrid.tri    = tri; 
%    disp('assigned tri in hgrid')    
end

ob.xy.flag  = 0;
ob.z.flag   = 0;
ob.xyz.flag = 0;
