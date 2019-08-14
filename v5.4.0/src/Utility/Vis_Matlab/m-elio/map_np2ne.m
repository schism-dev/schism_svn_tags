function [data_ne]=map_np2ne(gr,data_np)
%[data_ne]=map_np2ne(gr,data_np)
% gr -grid as returned by readGrid
% data_np data at nodes size(np,1)
% data_ne data at elem size(ne,1)
%
%Sergey Frolov march 2004

if size(data_np)~=gr.hgrid.np
  error('size(data_np)=[np,1]')
end

elem    = gr.hgrid.elem;

idx3    = find(elem(:,2)==3);
idx4    = find(elem(:,2)==4);

data_ne = zeros([gr.hgrid.ne 1]);

data_ne(idx3) = mean([data_np(elem(idx3,3)) data_np(elem(idx3,4)) data_np(elem(idx3,5))]')';
data_ne(idx4) = mean([data_np(elem(idx4,3)) data_np(elem(idx4,4)) data_np(elem(idx4,5)) data_np(elem(idx4,6))]')';
