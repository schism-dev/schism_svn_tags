function [grid]=gr_readGrid(grid,fname,grType)
%[grid]=gr_readGrid([grid,fname,grType=[hg|vg|sc|ec|dr]])
%generic intefcae to compile a grid structure from the files for
%horisontal|vertical|sideConnect grids
%
%Sergey Frolov march, 2004

if nargin == 1              %initialise
    [grid] = iniGrid(1);   
elseif nargin == 3
    if  grType == 'hg'    %horizontal grid
	grid.hgrid      = gr_readHGrid(fname);
        grid.hgrid.gtype =   'hg';
    elseif   grType == 'vg'    %vertical grid
	grid.vgrid      = gr_readVGrid(fname);
        grid.vgrid.gtype =   'vg';
    elseif   grType == 'sc'    %side centers connectivity table
	grid.sideCent   = readSCGrid(fname);
        grid            = gr_compute_lc(grid);
        grid.sideCent.gtype = 'sc';
    elseif   grType == 'ec'    %element centers
	grid.ecenters      = gr_readHGrid(fname);
        grid.ecenters.gtype = 'ec';        
    elseif   grType == 'dr'    %this is a pointer to a directory with sdtandart files
        grid = readGrDir(grid,fname);
    else
        error('you have to specify a grid of known type [hg|vg|sc|dr]')
    end
else
    error('wrong number of arguments')
end

% % %%%%%%%%%%%%%
% % %private functions
% % 
function [grid] = iniGrid(junk) 
    grid.hgrid.flag     = 0;
    grid.vgrid.flag     = 0;
    grid.sideCent.flag  = 0;
    grid.lc.flag        = 0;    

function grid = readGrDir(grid,fname)
    if exist(fname)~= 7
        error('with a "dr" option fname should be a directory')
    end
    fn  = fullfile(fname,'hgrid.gr3');
    if exist(fn) == 2
	grid.hgrid  = gr_readHGrid(fn);
        grid.hgrid.gtype =   'hg';
    else
        warning('hgrid.gr3 doesnt exist in the target directory')
    end
    fn  = fullfile(fname,'vgrid.in');
    if exist(fn) == 2
	grid.vgrid  = gr_readVGrid(fn);
        grid.vgrid.gtype =   'vg';        
    else
        warning('vgrid.in doesnt exist in the target directory')
    end
    fn  = fullfile(fname,'sidecenters_conn.bp');
    fn2 = fullfile(fname,'sidecenters.bp');
    fn3 = fullfile(fname,'sidecenters_snxy.bp');
    if exist(fn) == 2 & exist(fn2) == 2
        if exist(fn3) == 2 
	  grid.sideCent   = gr_readSCGrid(fn,fn2,fn3);
	else 
          grid.sideCent   = gr_readSCGrid(fn,fn2);
	end
        grid            = gr_compute_lc(grid);
        grid.sideCent.gtype = 'sc';        
    else
        warning('sidecenters_conn.bp or sidecenters.bp  doesnt exist in the target directory')        
    end
    fn  = fullfile(fname,'centers.bp');
    if exist(fn) == 2
	grid.ecenters  = gr_readHGrid(fn);
        grid.ecenters.gtype = 'ec';         
    else
        warning('centers.bp doesnt exist in the target directory')
    end
    
