function outVar = map_sb2z(gr,sb_d,flagLn,eta,flg99)
% outVar = map_sb2z(gr,sb_d,flagLn,[eta,flg99])
% interpolate selfe binary data to a single vertical Z level
% gr	- grid struct. or header of the selfe binary file
% sb_d 	- selfe binary data organized by levels [nLevs,np]
% flagLn- ['s'|'b'|number] 's' - surface layer 'b' - bottom layer
%		number	- layer at depth number 
% eta	- if class(flagLn)==double then also need eta - 
%		vector of elevation at each point [np,1]
% flg99	- flag to indicate what to do with -99 values in salinity and temperature
%	although i take care of most dry nodes ome of them may still 'survive'
%	if flag99==1 then -99 will be replaced by nan, otherwise ignored
%	default is flg99=0
% outVar- ouput vecotr [np,1]
%
% Sergey Frolov, May 2005

if nargin < 5
  flg99 = 0 ;
end

if size(sb_d)~=[gr.vgrid.nLevels gr.hgrid.np]
  error('size(sb_d) should be  [nLevs,np]')
end

if strcmp(class(flagLn),'char')&strcmp(flagLn,'s')
  outVar= sb_d(end,:)';
elseif strcmp(class(flagLn),'char')&strcmp(flagLn,'b')
  outVar= sb_d(1,:)';
elseif strcmp(class(flagLn),'double')
  zq	= double(flagLn);		%query depth
  dp 	= -gr.hgrid.depth;		%depth
  z	= sb_computeZlevels(-dp, eta, gr.vgrid);
        %index of feasible nodes
  fidx 	= find( (zq<=z(:,end))&(zq>=z(:,1))&z(:,end)~=dp );
  fz 	= z(fidx,:);			%feasible levels
  fd 	= sb_d(:,fidx)';		%feasible data
  ub 	= zeros(size(fz,1),1);		%will hold index of upper bounds
  fdi   = zeros(size(fz,1),1);		%feasible interpolated data
  for i = 2:size(z,2)
	%find upper bound braketing indx 
    ubIdx 	= find((fz(:,i)>zq)&(ub==0));
    ub(ubIdx) 	= i;
    fdi(ubIdx)= fd(ubIdx,i-1) + ...
        ( (fd(ubIdx,i)-fd(ubIdx,i-1))./(fz(ubIdx,i)-fz(ubIdx,i-1)) ).* ...
        (zq-fz(ubIdx,i-1));
  end
  outVar	= nan*zeros(size(z,1),1);
  outVar(fidx) 	= fdi;
else
  error('confused about flagLn')
end

if flg99
  idx         = find(outVar<=-98);
%  length(idx)
  outVar(idx) = nan;
end

