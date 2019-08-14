function [grid]=gr_compute_lc(grid)
%[grid]=gr_compute_lc(grid)
% compute coordinates of local axes
% grid - is a grid structure as returned by gr_readGrid + gr_compute_sc
%
%sergey frolov march, 2004
if grid.sideCent.flag
   %caution snx sny here are not the same as in Elcirc code
   %snxy here are just length of the triangles sides
	    %snx=x2-x1 (side vector x)
   snx = grid.hgrid.nodes(grid.sideCent.connect(:,3),2)-grid.hgrid.nodes(grid.sideCent.connect(:,2),2);
	    %sny=y2-y1 (side vector y)
   sny = grid.hgrid.nodes(grid.sideCent.connect(:,3),3)-grid.hgrid.nodes(grid.sideCent.connect(:,2),3);
        %normal vector as n=(-sny,snx)
   n   = -[-sny snx];   
        %tangential vector as y=(snx,sny)
   t   = [snx sny];   
        %normalise the vectors
   for i = 1:size(n,1)
       n(i,:) = n(i,:)/norm(n(i,:));
       t(i,:) = t(i,:)/norm(t(i,:));
   end
   grid.lc.n    = n;
   grid.lc.t    = t;
   grid.lc.flag = 1;
else
   error('sorry i nead side center connectivity information as in grid.sideCent.conn');
end
