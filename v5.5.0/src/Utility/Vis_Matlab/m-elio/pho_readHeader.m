function [ header ] = sz_readHeader( fname )
% [ header ] = sb_readHeader( fname )
% read header for the selfe binary file (fname)
% Data Format v5.00 (sz-levels)
%
% SF, Jan 2006, 
% SF March 2006, corrected bogus treatment of bIdx=0

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
header		    = computeStepIdx(header); 
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
	header.vgrid.nLevels   = fread(fid,1,'int32');
        header.vgrid.kz        = fread(fid,1,'int32');
	header.vgrid.h0	       = fread(fid,1,'float32');
        header.vgrid.hs        = fread(fid,1,'float32');
        header.vgrid.hc        = fread(fid,1,'float32');
        header.vgrid.theta_b   = fread(fid,1,'float32');
        header.vgrid.theta_f   = fread(fid,1,'float32');
        header.vgrid.zLevels   = fread(fid,header.vgrid.kz-1,'float32');
	header.vgrid.sLevels   = fread(fid,header.vgrid.nLevels-header.vgrid.kz+1,'float32');

function [header]=readHgrid(fid,header)
	header.hgrid.startPos   = ftell(fid);
	header.hgrid.type       = 'gr'; 		%it is a grid without a boundary
	header.hgrid.np         = fread(fid,1,'int32');
	header.hgrid.ne         = fread(fid,1,'int32');

        sp = ftell(fid);
	hgridTmp                = fread(fid,[4 header.hgrid.np],'float32')';
	header.hgrid.x          = hgridTmp(:,1);
	header.hgrid.y          = hgridTmp(:,2);
	header.hgrid.depth      = hgridTmp(:,3);

        fseek(fid,sp,'bof');
        hgridTmp                = fread(fid,[4 header.hgrid.np],'int32')';
        header.hgrid.bIdx       = hgridTmp(:,4);	%bottom level index

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
        bIdx = max(1,h.hgrid.bIdx);	%correct bIdx==0 for dry nodes
%        bIdx = h.hgrid.bIdx;	%with correction in readHgrid max(1,..) is absolete
        h.gridSize  = sum(h.vgrid.nLevels - bIdx+1);
    elseif h.flagDm == 2
        h.gridSize  = h.hgrid.np;
    else
        error('Confused about the dimension of the grid flagDm')
    end
    h.stepSize  = 2*4 + h.hgrid.np*4 + h.gridSize*4*h.flagSv;

function h = computeStepIdx(h)
% compute h.idx.idxNodes h.idx.idxLev
    if h.flagDm == 3 | h.flagDm == 2
        xv_x    = nan*ones([h.hgrid.np h.vgrid.nLevels]);
        xv_v    = nan*ones([h.hgrid.np h.vgrid.nLevels]);
        n_x     = [1:h.hgrid.np]';
        nvrt    = h.vgrid.nLevels;
        blx     = h.hgrid.bIdx;
        for i=1:nvrt
            idx           = find(blx<=i);
            xv_x(idx,i)   = n_x(idx,1);
            xv_v(idx,i)   = i*ones([length(idx) 1]);
        end
        xv_x            = xv_x';
        xv_v            = xv_v';
        idxNodes        = xv_x(:);
        h.idx.idxNodes  = uint32(idxNodes(find(~isnan(idxNodes))));
        idxLev          = xv_v(:);
        h.idx.idxLev    = uint8(idxLev(find(~isnan(idxLev))));
        h.idx.flag      = 1;
    else
        error('Confused about the dimension of the grid flagDm')
    end
% precompute mapping from selfe binary to selfe hotstart format
  h.idx.idx_all = zeros(h.hgrid.np,h.vgrid.nLevels);
  idx_eb=[1:length(h.idx.idxNodes)]';
  for i =1:h.vgrid.nLevels
    idxl        = find(h.idx.idxLev == i);
    idxn        = h.idx.idxNodes(idxl);
    h.idx.idx_all(idxn,i)  = idx_eb(idxl);
  end
  h.idx.idx_all=h.idx.idx_all(:);
  h.idx.idx_all_mask = logical(h.idx.idx_all~=0);

