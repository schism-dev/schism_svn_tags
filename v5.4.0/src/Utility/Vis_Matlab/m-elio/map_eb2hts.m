function [data_hts] = map_eb2hts(data_eb,flg)
%function [data_hts] = map_eb2hts(data_eb)
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

idxLev      = data_eb.h.idx.idxLev;
idxNodes    = data_eb.h.idx.idxNodes;
nlev        = data_eb.h.vgrid.nLevels;
np          = data_eb.h.hgrid.np;
data        = data_eb.data(:,1);

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
