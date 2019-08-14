function varOut = fort_readColumn(fname,n,prec)
%varOut = fort_readColumn(fname,n,prec)
%reads in the "fort" file where first two int32 are the dimensions, and the rest is a column major matrix of single precision
%read only column n 
% fname -filename
% n - column to read
% prec - output precision 'single' or 'double'
%
%sergey frolov august 2004

if nargin<3 
    prec='double'
end

fid         = fopen(fname,'r');
[mn,c]      = fread(fid,2,'int32');  %size m by n
fseek(fid, 8+mn(1)*(n-1)*4, 'bof');
if strcmp(prec,'single')
	[varOut,c]  = fread(fid,mn(1),'float32=>float32');
else
	[varOut,c]  = fread(fid,mn(1),'float32');
end
fid         = fclose(fid);