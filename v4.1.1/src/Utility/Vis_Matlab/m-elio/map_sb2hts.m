function [hts_data] = map_sb2hts(sb_h,sb_data,flg)
%function [hts_data] = map_sb2hts(sb_h,sb_data,flg)
% map a vector of data as extracted from slefe binary files dat(1:nLevs,1:np)
% to array as required by hotstart files dat([0|1]:nvrt,1:np). 
% this function is not as envovled as for elcirc, this one just reshpaes the array
%
% sb_h - header strucutre for selfe binary file
% sb_data - data read from selfe file
% hts_dat is a simple 2d aray
% flg = [0|1] 0-levels start from 0 index 1- levels start from 1, see elcirc for details
%
% Sergey Frolov May 2005

if nargin == 1
    flg =1;
end

if flg==1  %levels 1:nlev
  hts_data	= reshape(sb_data, sb_h.vgrid.nLevels, sb_h.hgrid.np);
else    %levels 0:nlev
  hts_data	= zeros(sb_h.vgrid.nLevels+1, sb_h.hgrid.np);
  hts_data(2:end,:) = reshape(sb_data, sb_h.vgrid.nLevels, sb_h.hgrid.np);
end
