function [data_sc]=map_np2sc(gr,data_np)
%[data_sc]=map_np2sc(grid,data_np)
% map data_np defined on grid nodes to data_sc defined on element side centers
% gr - is a grid structure (including gr.sideCent.connect)
% size(data_np) = [nlevs,gr.hgrid.np]; size(data_sc) = [nlevs,gr.sideCent.np];
%
% Sergey Frolov march 08, 2004
% vectorized on march 1, 2005

if ~gr.sideCent.flag
    error('sorry i need side center information as in grid.sideCent');
end
if size(data_np,2)~=gr.hgrid.np
  error('size(data_np,2) should be == to number of points in the grid')
end

data_sc = zeros([size(data_np,1) gr.sideCent.np]);
conn    = gr.sideCent.connect;

for i=1:size(data_np,1)
   data_sc(i,:)  = (data_np(i,conn(:,2))+data_np(i,conn(:,3)))/2;
end


