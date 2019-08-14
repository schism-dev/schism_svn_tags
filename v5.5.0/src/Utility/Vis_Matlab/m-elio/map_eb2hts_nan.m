function [data_hts] = map_eb2hts_nan(eb_h,eb_data)
%function [data_hts] = map_eb2hts_nan(eb_h, data_eb)
% map a vector of data as extracted from elcirc binary files dat(kbp:nvrt,1:np)
% to array as required by hotstart files dat([0|1]:nvrt,1:np). including filling in of the values above and below kfp kbp with nans
%
% eb_h is a header for the elcirc binary file
% eb_data is the data read by eb_readTimeStep
% data_hts is a simple 2d aray
%
% Sergey Frolov March 2004
% modified sf dec 2004

idxLev      = eb_h.idx.idxLev;
idxNodes    = eb_h.idx.idxNodes;
nlev        = eb_h.vgrid.nLevels;
np          = eb_h.hgrid.np;

data_hts = nan*zeros(np,nlev);
for i =1:nlev
    idxl        = find(idxLev == i);
    idxn        = idxNodes(idxl);
    data_hts(idxn,i)  = eb_data(idxl); 
end 
