function [ob]=ob_comp_wz(ob)	
% [ob]=ob_comp_wz(ob)
% computes weights and indexes for points ob.z to be interpolated from the grid ob.gr
% uses linear interpolation between bracating layers
%
% Sergey Frolov, April 2004

pz      = -ob.z.z;
pnp     = size(pz,1);
vgr     = (ob.gr.vgrid.zLevel - ob.gr.vgrid.zMsl);

for i=1:pnp
    ilb     = find(diff(vgr < pz(i)) == -1);   %lower bound index
    iub     = ilb +1;   %upper bound index
    w(i,2)  = (pz(i)-vgr(ilb))/(vgr(iub)-vgr(ilb));
    w(i,1)  = 1 - w(i,2);
    zi(i,:) = [ilb iub];
end

ob.z.pidx     = zi;
ob.z.w        = w;
ob.z.flag     = 2;
