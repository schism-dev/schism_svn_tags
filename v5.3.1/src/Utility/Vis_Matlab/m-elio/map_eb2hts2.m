function [data_hts] = map_eb2hts(h,data,flg)
%function [data_hts] = map_eb2hts2(h,data,flg)
% map a vector of data as extracted from elcirc binary files dat(kbp:nvrt,1:np)
% to array as required by hotstart files dat([0|1]:nvrt,1:np). including filling in of the values above and below kfp kbp
%
% data_eb is a structure as returned by [data_eb.h, data_eb.ts, data_eb.data]=eb_readTimeStep
% data_hts is a simple 2d aray
% flg = [0|1] 0-levels start from 0 index 1- levels start from 1, see elcirc for details
%
% Sergey Frolov March 2004

if nargin == 1
    flg =1;
end

idxLev      = h.idx.idxLev;
idxNodes    = h.idx.idxNodes;
nlev        = h.vgrid.nLevels;
np          = h.hgrid.np;

if flg  %levels 1:nlev
	data_hts = zeros(nlev, np);
    for i =1:nlev
        idxl        = find(idxLev == i);
        idxn        = idxNodes(idxl);
        data_hts(i,idxn)  = data(idxl)'; 
    end 
else    %levels 0:nlev
	data_hts = zeros(nlev+1, np);
    for i =1:nlev
        idxl        = find(idxLev == i);
        idxn        = idxNodes(idxl);
        data_hts(i+1,idxn)  = data(idxl)'; 
	end 
end
