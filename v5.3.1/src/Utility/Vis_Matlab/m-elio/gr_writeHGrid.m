function []=gr_writeHGrid(fname,grid)
%[status]=gr_writeHGrid(fname,grid)
%writes out a horisontal grid
%
%Sergey Frolov March, 2004

tic
if nargin <2 
    error('wrong number of inputs')
end
checkGrid(grid.hgrid) 

fid=fopen(fname,'w')
writeHeader(fid,grid.hgrid);
writeNodes(fid,grid.hgrid)
if grid.hgrid.type == 'gr'
   writeElem(fid,grid.hgrid);
elseif grid.hgrid.type == 'gb'
    writeElem(fid,grid.hgrid);
    writeEofLines(fid,grid.hgrid);
end
fclose(fid);
toc

% disp(['read grid file ' grid.fname ' in ' num2str(toc) ' seconds'])


%%%%%%%%%%%%%
%private functions

function [] = checkGrid(grid) 
	if grid.type=='bp'
        if [grid.np 4] ~= size(grid.nodes), error('size(nodes) ~= [np 4]');end
        disp('writing a build point file')
    elseif grid.type=='gr'|'gb'
        if [grid.np 4] ~= size(grid.nodes), error('size(nodes) ~= [np 4]');end
        if [grid.ne 6] ~= size(grid.elem), error('size(elem) ~= [ne 6]');end
        disp('writing a grid file')
    else
        disp('unkown grid type')
    end

function []=writeHeader(fid,grid)
	fprintf(fid,'%s \n',grid.name);
    if grid.type == 'bp' 
        lineOut = num2str([grid.np]);
	elseif grid.type == 'gb'|'gr' 
        lineOut = num2str([grid.ne grid.np]);        
	else
        error('oops, what kind of file is this?')
	end
	fprintf(fid,'%s \n',lineOut);
    
function []=writeNodes(fid,grid)
    fprintf(fid,'%u %6.6f %6.6f %g \n',grid.nodes');
    
function []=writeElem(fid,grid)
    fprintf(fid,'%u %u %u %u %u %u \n',grid.elem');
    
function []=writeEofLines(fid,grid)    
    for i=1:length(grid.eofLines)
		fprintf(fid,'%s \n',grid.eofLines{i});
    end
