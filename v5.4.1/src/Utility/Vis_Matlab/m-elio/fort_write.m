function status = fort_write(fname, varIn)
%status = fort_write(fname, varIn)
%writes a "fort" file
%see code for file format descr
%
%sergey frolov august 2004

fid     = fopen(fname,'w');
c       = fwrite(fid,size(varIn),'int32');
c       = fwrite(fid,varIn,'float32');
fclose(fid);
status =1;
