function fid = sb_readHeader( fname, h, fidFlag )
% [fid] = sb_readHeader( fname, h, [fidFlag=0] )
% write header (h) to the selfe binary file (fname)
% if fidFalg return open fid
% fid - fid of the file 0 if fidFlag=0
%
% SF, November 2005, 

% tic
fid = fopen(fname, 'w');
if nargin <3
  fidFlag = 0;
end

%hheader
c = fwrite(fid,h.dataFormat,'uchar');
c = fwrite(fid,h.version,'uchar');
c = fwrite(fid,h.startTime,'uchar');
c = fwrite(fid,h.varType,'uchar');
c = fwrite(fid,h.varDimension,'uchar');
c = fwrite(fid,h.nSteps,'int32');
c = fwrite(fid,h.dt,'float32');
c = fwrite(fid,h.skip,'int32');
c = fwrite(fid,h.flagSv,'int32');
c = fwrite(fid,h.flagDm,'int32');

%vgrid
c = fwrite(fid,h.vgrid.ivcor,'int32');
c = fwrite(fid,h.vgrid.h0,'float32');
c = fwrite(fid,h.vgrid.hc,'float32');
c = fwrite(fid,h.vgrid.theta_b,'float32');
c = fwrite(fid,h.vgrid.theta_f,'float32');
c = fwrite(fid,h.vgrid.nLevels,'int32');
c = fwrite(fid,h.vgrid.sLevels,'float32');

%hgrid
c = fwrite(fid,h.hgrid.np,'int32');
c = fwrite(fid,h.hgrid.ne,'int32');
hgridTmp = [h.hgrid.x h.hgrid.y h.hgrid.depth];
c = fwrite(fid,hgridTmp','float32')';
hgridTmp = h.hgrid.elem(:,2:5);
c = fwrite(fid,hgridTmp','int32')';

%toc
if ~fidFlag
  fid = fclose(fid);
end
