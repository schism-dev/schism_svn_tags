function [z] = sb_computeZlevels(dp, eta, vgrid)
%[zLevels] = sb_computeZlevels(dp, eta, vgrid)
% convert S, sigma, mixed coordinates to Z levels
% dp 	- depth 	[np,1]
% eta 	- elevation	[np,1]
%	np is number of points, 
%	points don't need to coinside with points of the horizontal grid
% vgrid - vertical grid structurte (e.g. as returned in sb_header)
% zLevel- [np,nLevs] Z levle at each  
%	zLevs(:,1)=dp, zLevs(:,end)=eta
%
% NOTE: there may be a suttle difference in the output of these function 
%	and zLevels in ?_zcor.63. 
%	for dry nodes ?_zcor.63 has zLevs=0, I have zLevs=dp
%	hence you can find later almost all dry nodes as idx=find(zLevs(:,end)==dp)
%
% Sergey Frolov May 2005

z       = zeros(size(eta,1),length(vgrid.sLevels));
for i=1:length(vgrid.sLevels)
        %compute C(sigma)
  sigma         = vgrid.sLevels(i);
  C             = (1-vgrid.theta_b)*sinh(vgrid.theta_f*sigma)/sinh(vgrid.theta_f) + ...
        vgrid.theta_b*( tanh(vgrid.theta_f*(sigma+0.5))-tanh(vgrid.theta_f/2) )/ ...
                (2*tanh(vgrid.theta_f/2));
        %special case where dp<=hc
  idx           = (dp<=vgrid.hc);
  z(idx,i)      = sigma*(dp(idx)+eta(idx))+eta(idx);
        %normal case
  z(~idx,i)     = eta(~idx)*(1+sigma) + vgrid.hc*sigma + (dp(~idx)-vgrid.hc)*C;
        %special case 2
  criterion2    = -vgrid.hc-(dp-vgrid.hc)*vgrid.theta_f/sinh(vgrid.theta_f);
  idx           = find(eta<=criterion2);
  if length(idx)>0
%    warning('S->Z coordinate transformation may not be valid')
    etam        = 0.98*( -vgrid.hc-(dp-vgrid.hc)*vgrid.theta_f/sinh(vgrid.theta_f) );
    sigma_hat   = (z(idx,i)-etam(idx))./(dp(idx)+etam(idx));
    z(idx,i)    = sigma_hat.*(dp(idx)+eta(idx))+eta(idx);
  end
        %special case of dry nodes
        %zoseph assigns z=0 in this case,
        %i think it's more consistent to assign depth value in this case
  idx           = find((dp + eta)<=vgrid.h0);
%  z(idx,i)      = 0;           %zhoseph
  z(idx,i)      = dp(idx);      %sergey

end

