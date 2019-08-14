function varOut = fort_readColumns(fname,nn,prec)
%varOut = fort_readColumns(fname,nn,prec)
%reads in the "fort" file where first two int32 are the dimensions, and the rest is a column major matrix of single precision
%read only column n 
% fname -filename
% nn - vector of column numbers to read
% prec - output precision 'single' or 'double'
%
%sergey frolov august 2004

if nargin<3 
    prec='double';
end

fid         = fopen(fname,'r');
[mn,c]      = fread(fid,2,'int32');  %size m by n
if max(nn)>mn(2)
    error(['column number ' num2str(max(nn)) ' exceeds number of columns ' num2str(mn(2)) ' in file ' fname]);
end

if strcmp(prec,'single')
	varOut      = single(zeros(mn(1), length(nn)));
else
	varOut      = zeros(mn(1), length(nn));
end

for i = 1:length(nn)
    n=nn(i);
	fseek(fid, 8+mn(1)*(nn(i)-1)*4, 'bof');
	if strcmp(prec,'single')
		[varOut(:,i),c]  = fread(fid,mn(1),'float32=>float32');
	else
		[varOut(:,i),c]  = fread(fid,mn(1),'float32');
	end
end

fid         = fclose(fid);
