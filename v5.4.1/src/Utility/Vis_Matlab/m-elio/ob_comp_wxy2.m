function [ob,tidx]=ob_comp_wxy2(ob)	
% [ob,tidx]=ob_comp_wxy2(ob)
% computes weights and indexes for points ob.xy to be interpolated from the grid ob.gr
%handles case when some points are on land
%
% Sergey Frolov April 2004
%modified Sergey Frolov December 2004
%keyboard

px      = ob.xy.x;
py      = ob.xy.y;
pnp     = size(px,1);
hgr     = ob.gr.hgrid;
tidx    = tsearch(hgr.x,hgr.y,hgr.tri,px,py);
nanidx	= isnan(tidx);
if sum(nanidx)
%    error('cant proceed, few of your points are outside of the grid')
  pidx  = hgr.tri(tidx(~nanidx),:);
  px	= px(~nanidx);
  py    = py(~nanidx);
  ob.xy.x = ob.xy.x(~nanidx);
  ob.xy.y = ob.xy.y(~nanidx);
  ob.xy.tridx = tidx(~nanidx);
else
  pidx    = hgr.tri(tidx,:); 
  ob.xy.tridx = tidx;
end
%!!! important order points conterclockwise !!!
pidx    = pidx(:,[3 1 2]);  

tx  = hgr.x(pidx);
ty  = hgr.y(pidx);
if pnp==1  %rotate dimensions
    tx = tx';
    ty = ty';
end

%this section of code is borrowed from ELIO v1.33
aum = tx(:,2).* ty(:,3) - tx(:,3).* ty(:,2);
bum = ty(:,2) - ty(:,3);
cum = tx(:,3) - tx(:,2);
ado = tx(:,3).* ty(:,1) - tx(:,1).* ty(:,3);
bdo = ty(:,3) - ty(:,1);
cdo = tx(:,1) - tx(:,3);
atr = tx(:,1).* ty(:,2) - tx(:,2).* ty(:,1);
btr = ty(:,1) - ty(:,2);
ctr = tx(:,2) - tx(:,1);
arei = (aum + ado + atr);
w(:,3) = (atr + btr.* px + ctr.* py)./arei;
w(:,2) = (ado + bdo.* px + cdo.* py)./arei;
w(:,1) = 1.0 - w(:,2) - w(:,3);


%assignment of vars in ob
%ob.xy.tridx     = tidx;
ob.xy.pidx      = pidx;
ob.xy.w         = w;
ob.xy.flag      = 2; 
