function [dataSc]=map_np2sc(grid,dataNp)
%[dataSc]=map_elem2sc(grid,dataEl)
% map dataNp defined on grid nodes to dataSc on element side centers
% grid - is a grid structure as returned by readGrid + computeSc
%
% Sergey Frolov march 08, 2004

if grid.sideCent.flag
	conn    = grid.sideCent.connect;
	dataSc  = (dataNp(conn(:,2))+dataNp(conn(:,3)))/2;
else
    error('sorry i need side center information as in grid.sideCent');
end