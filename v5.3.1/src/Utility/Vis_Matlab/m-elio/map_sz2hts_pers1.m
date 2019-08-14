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
% fast version with persistent vars, SF, Apr 2006

%using persistent varaible to precompute the loop
persistent idxl_all idxn_all 
if isempty(idxl_all)
  for i =1:eb_h.vgrid.nLevels
    idxl_all{i} = (eb_h.idx.idxLev == i);
    idxn_all{i} = eb_h.idx.idxNodes(idxl_all{i});
  end
%else
  %do some smart checking if we are working with the same hgrid/vgrid combination
end

%output format with nLevs = nvrt if flg =1; else flg=0 -> nLevs = nvrt+1
if nargin <3
  flg=1;
end
if flg == 1
  data_hts = nan*zeros(eb_h.hgrid.np,eb_h.vgrid.nLevels);
  iplus=0;
elseif flg==0
  data_hts = nan*zeros(eb_h.hgrid.np,eb_h.vgrid.nLevels+1);
  iplus=1;
else
  error('unknown flg')
end

%core re-assignment loop
for i =1:eb_h.vgrid.nLevels
  ii = i+iplus;
  data_hts(idxn_all{i},ii)  = eb_data(idxl_all{i}); 
end 


