function [vn,vt]=map_gc2lc(gr,ugsc, vgsc)
%[vn,vt]=map_gc2lc(gr,ugsc, vgsc)
% converts data defined at global coordinates (ugsc, vgsc) to local coordinates (vn, vt)
% size(ugsc)==size(vgsc)==[num_levels,num_sides]
% gr - a comp grid as returned by readGrid + computeSc +computeLocalCoord
%
% sergey frolov march 2005

if size(ugsc)~=size(vgsc)
  error('size vn and vt should be the same')
elseif size(vgsc,2)~=gr.sideCent.np
  error('size(vn,2) should be == to number of sides in the grid')
end

vn = zeros(size(vgsc));
vt = zeros(size(vgsc));

snx = gr.sideCent.snx';
sny = gr.sideCent.sny';

for i=1:size(vn,1)
    vn(i,:)=  ugsc(i,:).*snx+vgsc(i,:).*sny;
    vt(i,:)= -ugsc(i,:).*sny+vgsc(i,:).*snx;
end

