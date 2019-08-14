function [ugsc, vgsc]=map_lc2gc(gr,vn,vt)
%[ugsc, vgsc]=map_lc2gc(gr,vn,vt)
% converts data defined at local coordinates (vn, vt) to global coordinates (ugsc, vgsc)
% size(vn)==size(vt)==[num_levels,num_sides]
% gr - a comp grid as returned by readGrid + computeSc +computeLocalCoord
%
%sergey frolov march 2005


if size(vn)~=size(vt)
  error('size vn and vt should be the same')
elseif size(vn,2)~=gr.sideCent.np
  error('size(vn,2) should be == to number of sides in the grid')
end

ugsc = zeros(size(vn));
vgsc = zeros(size(vn));

snx = gr.sideCent.snx';
sny = gr.sideCent.sny';

for i=1:size(vn,1)
    ugsc(i,:) = vn(i,:).*snx-vt(i,:).*sny;
    vgsc(i,:) = vn(i,:).*sny+vt(i,:).*snx;
end

