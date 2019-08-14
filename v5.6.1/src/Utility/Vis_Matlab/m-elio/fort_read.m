function varOut = fort_read(fname,prec)
%varOut = fort_read(fname,[prec])
%reads in the "fort" file where first two int32 are the dimensions, and the rest is a column major matrix of single precision
% fname -filename
% prec - output precision 'single' or 'double'
%
%sergey frolov august 2004

if nargin<2 
    prec='double'
end

fid         = fopen(fname,'r');
[mn,c]      = fread(fid,2,'int32');  %size m by n
if strcmp(prec,'single')
	[varOut,c]  = fread(fid,mn','float32=>float32');
else
	[varOut,c]  = fread(fid,mn','float32');
end
fid         = fclose(fid);
