function [ob]= ob_loadSta(ob, xyzFname)
%[ob]= ob_loadSta(ob, xyzFname)
% loads trasect file (*.bp) into ob datastructure
%INPUT/OUTPUT
% ob - observation datastructure as initialized by ob_ini
% xyzFanme - xyz file of the format: ('%8c %f %f %f') = (name x y z)
%
% Sergey Frolov, Dec 2004

if exist(xyzFname)~=2
    error(['file ' xyzFname ' doesnt exist'])
end

%load points
[ob.xyz.name,ob.xy.x,ob.xy.y,ob.z.z]=textread(xyzFname, '%8c %f %f %f');
ob.xy.flag  = 1;    %point loaded
ob.z.flag   = 1;    %point loaded
ob.xyz.flag = 1;    %point loaded
display('xyz points loaded')

%%%%
%xy
%compute interpolation weights
[ob]        = ob_comp_wxy(ob);
ob.xy.flag  = 2;    %weights computed
display('xy weights computed')
%assemble sparce observation matrix H
[ob]=ob_obxy2h(ob);
ob.xy.flag  = 3;    %matrix H assembled
display('xy matrix H assembled')

%%%%
%z
%compute interpolation weights
[ob]        = ob_comp_wz(ob);
ob.z.flag   = 2;    %weights computed
display('z weights computed')
%assemble sparce observation matrix H
[ob]        = ob_obz2h(ob);
ob.z.flag   = 3;    %matrix H assembled
display('z matrix H assembled')
