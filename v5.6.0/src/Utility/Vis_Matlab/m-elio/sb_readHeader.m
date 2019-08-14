function [ header ] = sb_readHeader( fname )
% [ header ] = sb_readHeader( fname )
% read header for the selfe binary file (fname)
%
% SF, May 2005, 

% tic
header.fname        = fname;
header.stepSize     = 0;

fid                 = fopen(fname);
header.fid          = fid;
header              = readHHeader(fid,header);
header              = readVgrid(fid,header);
header              = readHgrid(fid,header);
header.dataStartPos = ftell(fid);
header 		    = computeStepSize(header);	 
header.fid          = fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%private functions
function [header]=readHHeader(fid,header)
	header.dataFormat   = char(fread(fid,48,'uchar')');
	header.version      = char(fread(fid,48,'uchar')');
	header.startTime    = char(fread(fid,48,'uchar')');
	header.varType      = char(fread(fid,48,'uchar')');
	header.varDimension = char(fread(fid,48,'uchar')');
	header.nSteps       = fread(fid,1,'int32');
	header.dt           = fread(fid,1,'float32');
	header.skip         = fread(fid,1,'int32');
	header.flagSv       = fread(fid,1,'int32');
	header.flagDm       = fread(fid,1,'int32');
%sb	header.zDes         = fread(fid,1,'float32');

function [header]=readVgrid(fid,header)
	header.vgrid.startPos  = ftell(fid);
	header.vgrid.ivcor     = fread(fid,1,'int32');
	header.vgrid.h0	       = fread(fid,1,'float32');
        header.vgrid.hc        = fread(fid,1,'float32');
        header.vgrid.theta_b   = fread(fid,1,'float32');
        header.vgrid.theta_f   = fread(fid,1,'float32');
        header.vgrid.nLevels   = fread(fid,1,'int32');
	header.vgrid.sLevels   = fread(fid,header.vgrid.nLevels,'float32');

function [header]=readHgrid(fid,header)
	header.hgrid.startPos   = ftell(fid);
	header.hgrid.type       = 'gr'; %it is a grid without a boundary
	header.hgrid.np         = fread(fid,1,'int32');
	header.hgrid.ne         = fread(fid,1,'int32');
	hgridTmp                = fread(fid,[3 header.hgrid.np],'float32')';
	header.hgrid.x          = hgridTmp(:,1);
	header.hgrid.y          = hgridTmp(:,2);
	header.hgrid.depth      = hgridTmp(:,3);
	header.hgrid.nodes	= [[1:header.hgrid.np]' header.hgrid.x header.hgrid.y header.hgrid.depth];
        header			= readElem_3(fid,header);

function [header]=readElem_3(fid,header)
    ne                  = header.hgrid.ne;
    elem                = nan*ones(ne,5);
    for i=1:ne
        toq = fread(fid,1,'int32');
        elem(i,1)       = toq;
        if toq == 3
            elem(i,2:4) = fread(fid,3,'int32')';
        elseif toq == 4
            elem(i,2:5) = fread(fid,4,'int32')';
        else
            error(['oops problems with reading element ' num2str(i)])
        end
    end
    header.hgrid.elem   = [[1:size(elem,1)]' elem];;
    
function h = computeStepSize(h)
    if h.flagDm == 3
        h.gridSize  = h.hgrid.np*h.vgrid.nLevels;
    elseif h.flagDm == 2
        h.gridSize  = h.hgrid.np;
    else
        error('Confused about the dimension of the grid flagDm')
    end
    h.stepSize  = 2*4 + h.hgrid.np*4 + h.gridSize*4*h.flagSv;

