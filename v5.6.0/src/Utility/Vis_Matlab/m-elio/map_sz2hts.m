function [data_hts] = map_sz2hts(eb_h,eb_data,flg)
%function [data_hts] = map_sz2hts(eb_h, data_eb)
% map a vector of data as extracted from elcirc or selfe (sz) binary files dat(kbp:nvrt,1:np)
% to array as required by hotstart files dat([0|1]:nvrt,1:np). including filling in of the values above and below kfp kbp with nans
%
% eb_h is a header for the binary file
% eb_data is the data read by ??_readTimeStep
% data_hts is a simple 2d aray
%
% Sergey Frolov March 2004
% modified from map_eb2hts_nan by SF, Jan 2006


data_hts = nan(size(eb_h.idx.idx_all));
data_hts(eb_h.idx.idx_all_mask)=eb_data(eb_h.idx.idx_all(eb_h.idx.idx_all_mask));
data_hts = reshape(data_hts,eb_h.hgrid.np,eb_h.vgrid.nLevels);
%data_hts(~eb_h.idx.idx_all_mask)=nan;
