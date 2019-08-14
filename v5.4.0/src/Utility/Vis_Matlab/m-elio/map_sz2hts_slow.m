function [data_hts] = map_sz2hts(eb_h,eb_data,flg)
%function [data_hts] = map_sz2hts(eb_h, data_eb,[flg=1])
% map a vector of data as extracted from elcirc or selfe (sz) binary files dat(kbp:nvrt,1:np)
% to array as required by hotstart files dat([0|1]:nvrt,1:np). including filling in of the values above and below kfp kbp with nans
%
% eb_h is a header for the binary file
% eb_data is the data read by ??_readTimeStep
% data_hts is a simple 2d aray
%
% Sergey Frolov March 2004
% modified from map_eb2hts_nan by SF, Jan 2006

if nargin <3
  flg=1;
end

idxLev      = eb_h.idx.idxLev;
idxNodes    = eb_h.idx.idxNodes;
nlev        = eb_h.vgrid.nLevels;
np          = eb_h.hgrid.np;

if flg == 1
  data_hts = nan*zeros(np,nlev);
elseif flg==0
  data_hts = nan*zeros(np,nlev+1);
else
  error('unknown flg')
end

for i =1:nlev
  idxl        = find(idxLev == i);
  idxn        = idxNodes(idxl);
  if flg == 1
    ii=i;
  elseif flg==0
    ii=i+1;
  end
  data_hts(idxn,ii)  = eb_data(idxl); 
end 


