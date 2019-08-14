function mn = fort_readHEad(fname)
%mn = fort_readHEad(fname)
%reads the header of the fort file
% fname -filename
% mn - [numRows numColumns]
%
%sergey frolov august 2004


fid         = fopen(fname,'r');
[mn,c]      = fread(fid,2,'int32');  %size m by n
fid         = fclose(fid);