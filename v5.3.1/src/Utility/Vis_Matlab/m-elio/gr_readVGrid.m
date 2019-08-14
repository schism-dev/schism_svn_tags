function [grid]=gr_readVGrid(fname)
%[grid]=gr_readVGrid(fname)
%reads in the vertical grid in to a grid structure

fid     = fopen(fname);

line    = sscanf(fgetl(fid),'%i %g',2);
nLevels = line(1);
zMsl    = line(2);
levels  = fscanf(fid,'%i %g %g',[3 nLevels ])';

grid.nLevels    = nLevels;
grid.zMsl       = zMsl;
grid.levels     = levels;
grid.zLevel     = levels(:,3);

grid.flag       = 1;

fid=fclose(fid);
