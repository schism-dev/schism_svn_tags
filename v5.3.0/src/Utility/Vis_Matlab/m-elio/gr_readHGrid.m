function [grid]=gr_readHGrid(fname)
%[grid]=gr_readHGrid(fname)
%reads in the gr3 grid or build point files in to a grid structure
%
% Sergey Frolov March, 2004

tic

lines   = textread(fname,'%s','delimiter','\n','whitespace','');
lines   = strrep(lines,'D+00','e+00');        %fortran can output expenet as D sometines ????

grid    = iniGrid(fname);
grid    = readHeader(lines,grid);
grid    = readNodes(lines,grid);
if grid.type == 'gb'
    grid    = readElem(lines,grid);
    grid    = readEofLines(lines,grid);
else
    grid    = readEofLines(lines,grid);        
end
grid.flag = 1;

%disp(['read grid file ' grid.fname ' in ' num2str(toc) ' seconds'])


%%%%%%%%%%%%%
%private functions

function [grid] = iniGrid(fname) 
	grid.fname  = fname;
	grid.name   = [];
	grid.np     = [];
	grid.ne     = [];
	grid.nodes  = [];
	grid.elem   = [];
	grid.eofLines       = {};
    grid.nn=[];grid.x=[];grid.y=[];grid.depth=[];

function [grid]=readHeader(lines,grid)
	grid.name   = lines{1};
	tmp         = sscanf(lines{2},'%i');
	if length(tmp) == 1
        np      = tmp(1);
        ne      = 0;
        type    = 'bp';     %build points
        disp('this is a build point file')
	elseif length(tmp) == 2
        ne      = tmp(1);
        np      = tmp(2);
        type    = 'gb';     %grid with possibly a boundary info
        disp('this is a grid file')
	else
        error('error reading ne np')
	end
	grid.type = type;
    grid.ne   = ne;
    grid.np   = np;

function [grid]=readNodes(lines,grid)
    np      = grid.np;
    nodes   = nan*ones(np,4);
    for i=1:np
        nodes(i,:) = sscanf(lines{i+2},'%i %f %f %f')';
    end
    grid.nodes = nodes;
    grid.nn=nodes(:,1); grid.x=nodes(:,2); grid.y=nodes(:,3); grid.depth=nodes(:,4);
    
function [grid]=readElem(lines,grid)
    ne      = grid.ne;
    np      = grid.np;
    elem    = nan*ones(ne,6);
	for i=1:ne
        line     = sscanf(lines{np+2+i},'%i %i %i %i %i %i')';
        if length(line) == 6
            elem(i,1:6) = line;
        elseif length(line) == 5
            elem(i,1:5) = line;
        else
            warning(['error reading element ' num2str(i)]);
        end
	end
	grid.elem   = elem;

function [grid]=readEofLines(lines,grid)  
    bndStart = 2+grid.np+grid.ne+1;
    if length(lines) > bndStart
		grid.eofLines = lines(2+grid.np+grid.ne+1:end);
    else
        grid.eofLines   = {};
        grid.type       = 'gr';     %grid with no boundary info
    end
