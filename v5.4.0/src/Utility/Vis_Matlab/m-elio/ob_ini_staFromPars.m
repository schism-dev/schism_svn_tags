function [ob itgFlag]= ob_ini_staFromPars(gr, xyz, xyzType)
%[ob]= ob_ini_staFromPars(gr, xyz, xyzType)
% loads trasect file (*.bp) into ob datastructure
%INPUT/OUTPUT
% gr 	- grid strucutre with hgrid and vgrid fileds
% xyz	- xyz vector [x y z] 
% xyzType- ['xy','xyz'] type of the observation
%
% Sergey Frolov, Jun 2005

%initialise ob strucutre
ob      = ob_ini(gr.hgrid, gr.vgrid);

%assign xy points
ob.xy.x     = xyz(:,1);
ob.xy.y     = xyz(:,2);
ob.xy.flag  = 1;    %points loaded
[ob,itgFlag]        = ob_comp_wxy(ob);
ob.xy.flag  = 2;    %weights computed
[ob]=ob_obxy2h(ob);
ob.xy.flag  = 3;    %matrix H assembled
%display('xy matrix H assembled')

%%%%
%z
if strcmp(xyzType, 'xyz')
  ob.z.z      = xyz(3);
  [ob]        = ob_comp_wz(ob);
  ob.z.flag   = 2;    %weights computed
  display('z weights computed')
  %assemble sparce observation matrix H
  [ob]        = ob_obz2h(ob);
  ob.z.flag   = 3;    %matrix H assembled
%  display('z matrix H assembled')
end
