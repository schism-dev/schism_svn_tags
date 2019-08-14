function [gr] =gr_tri1(gr)
% [gr] =gr_tri1(gr)
% triangulates grids in a structure gr (as returned by readGrid)
%only hgrid
% 
% Sergey Frolov April 2004

if gr.hgrid.flag ==1 & gr.hgrid.ne ==0
    grid        = gr.hgrid;
    tri         = delaunay(grid.nodes(:,2),grid.nodes(:,3));
    grid.ne     = size(tri,1);
    grid.elem   = [[1:grid.ne]', 3*ones(grid.ne,1), tri, nan*ones(grid.ne,1)];
    gr.hgrid    = grid;
%    disp('triangulated hgrid')
end

if gr.hgrid.flag ==1 
    grid        = gr.hgrid;
    tri         = delaunay(grid.nodes(:,2),grid.nodes(:,3));
    grid.tri    = tri; 
    gr.hgrid    = grid;   
%    disp('assigned tri in hgrid')    
end
