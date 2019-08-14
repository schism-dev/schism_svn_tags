function [h,ts,data] = eb_readTimeStep(h,n)
%[h,ts,data] = eb_readTimeStep(h,n)
% read one time step (n) from elcirc binary files 
% h -header structure as returned by readHeader
% n - array of timestep numbers 
%
% sergey frolov march 8, 2004

tic
if n>h.nSteps
    error('time step out of range')
end
if ~h.fid
    fid     = fopen(h.fname);
    h.fid   = fid;
end
if ~h.stepSize
    h = computeStepSize(h);
end
if ~h.idx.flag
	h = computeStepIdx(h);
end

data=zeros([h.gridSize h.flagSv length(n)]);
for i=1:length(n)
% 	[ts{i},d(:,:,i)]=readTs(h,n(i));
 		[t,d]=readTs(h,n(i));
        data(:,:,i)=d;
        ts{i}=t;
end

h.fid=fclose(h.fid);
% disp(['Read ' num2str(length(n)) ' time steps in ' num2str(toc)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%private functions
function h = computeStepSize(h)
	if h.flagDm == 3
        h.gridSize  = sum(h.vgrid.nLevels - h.hgrid.bottomLayer+1);
    elseif h.flagDm == 2
        h.gridSize  = h.hgrid.np;
    else
        error('Confused about the dimension of the grid flagDm')
    end
    h.stepSize  = 2*4 + h.hgrid.np*4 + h.gridSize*4*h.flagSv;
    
function [ts,data]=readTs(h,n)
    fseek(h.fid, h.dataStartPos + (h.stepSize)*(n-1), -1);
    
    ts.t        = fread(h.fid,1,'float32');
    ts.tit      = fread(h.fid,1,'int32');
	ts.sidx     = uint8(fread(h.fid, h.hgrid.np, 'int32'));
%     ts          = computeFlagSurf(h,ts);
%     disp(['timestep ' num2str(ts.tit) ' at time ' num2str(ts.t)])
    
	data    = fread(h.fid, [h.flagSv h.gridSize], 'float32')';            
    
function h = computeStepIdx(h)
    if h.flagDm == 3
		xv_x    = nan*ones([h.hgrid.np h.vgrid.nLevels]);
		xv_v    = nan*ones([h.hgrid.np h.vgrid.nLevels]);
		n_x     = [1:h.hgrid.np]';
		nvrt    = h.vgrid.nLevels;      
		blx     = h.hgrid.bottomLayer;
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
	elseif h.flagDm == 2
		h.idx.idxNodes  = uint32([1:1:h.hgrid.np]');
        h.idx.idxLev    = [];
        h.idx.flag      = 1;
    else
        error('Confused about the dimension of the grid flagDm')
    end

%obsolleet, can have better implentation with h.idx.idxNodes
function ts = computeFlagSurf(h,ts)
    if h.flagDm == 3
        flagSurf = nan*ones([h.hgrid.np h.vgrid.nLevels]);
		blx      = h.hgrid.bottomLayer;    
		sidx     = ts.sidx;
        nvrt     = h.vgrid.nLevels;
		for i=1:nvrt
            idx             = find(blx<=i);
            flagSurf(idx,i) = ones([length(idx) 1]);
            idx             = find(sidx<i);
            flagSurf(idx,i) = zeros([length(idx) 1]);
		end
        flagSurf    = flagSurf';
        flagSurf    = flagSurf(:);
        ts.flagSurf = flagSurf(find(~isnan(flagSurf)));
    elseif h.flagDm == 2
        ts.flagSurf = ones([h.gridSize 1]);
    else
        error('Confused about the dimension of the grid flagDm')
    end