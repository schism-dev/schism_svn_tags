function [data_lc]=map_gc2lc(hg,h,data_gc)
%[data_lc]=gc2lc(hg,data_gc)
% converts data defined at global coordinates to local coordinates
% data_gc[x,y] data_lc[n,t]
% hg - a comp grid as returned by readGrid + computeSc +computeLocalCoord
%
%sergey frolov march 2004

tic
n       = hg.lc.n;
t       = hg.lc.t;
data_lc = zeros(size(data_gc));
idx     = h.idx.idxNodes;

for i=1:size(data_gc,1)
    data_lc(i,1) = data_gc(i,:)*n(idx(i),:)';
    data_lc(i,2) = data_gc(i,:)*t(idx(i),:)';    
end

toc
