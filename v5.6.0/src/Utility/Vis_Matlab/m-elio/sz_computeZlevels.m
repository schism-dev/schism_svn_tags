function [z] = sz_computeZlevels(dp, eta, vgrid)
%[zLevels] = sz_computeZlevels(dp, eta, vgrid)
% convert sigma-z coordinates to Z levels,
% also works for pure sigma coordinates
% dp 	- depth 	[np,1]
% eta 	- elevation	[np,1]
%	np is number of points, 
%	points don't need to coinside with points of the horizontal grid
% vgrid - vertical grid structurte (e.g. as returned in sb_header)
% zLevel- [np,nLevs] Z level at each  
%	zLevs(:,1)=dp, zLevs(:,end)=eta
%
% NOTE: there may be a suttle difference in the output of this function 
%	and zLevels in ?_zcor.63. 
%	for dry nodes ?_zcor.63 has zLevs=0, I have zLevs=dp
%	hence you can find later almost all dry nodes as idx=find(zLevs(:,end)==dp)
%
% Sergey Frolov May 2005


if isfield(vgrid,'zLevels')	%sigma-z

  %1) compute sigma levels reusing code for the pure sigma levles
  %   includes faking dp=vgrid.hs if dp>hs

  [zs] = sb_computeZlevels(min(vgrid.hs,dp), eta, vgrid);

  %2) array of pure z levels
  zz = repmat(vgrid.zLevels',size(zs,1),1);
  
  %3) set the depth of the bottom layer to bathymetry
  dp2=-repmat(dp,1,size(zz,2)+1);
  tmp1=[zz -vgrid.hs*ones(size(dp))]<=dp2;
  tmp2=tmp1(:,1:end-1)-tmp1(:,2:end);
  zz(tmp2==1)=dp2(tmp2==1);
  
  %4)concatinate zz and zs layers
  z=[zz zs];

else	%pure sigma
  z = sb_computeZlevels(dp, eta, vgrid);
end
