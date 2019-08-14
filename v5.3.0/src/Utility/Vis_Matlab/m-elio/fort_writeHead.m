function status = fort_writeHEad(fname, mn)
%status = fort_writeHead(fname, mn)
%writes a header (size of the array) into a "fort" file
%see code for file formate descr
%
%sergey frolov august 2004

fid     = fopen(fname,'w');
if fid ==0 
    error(['failed openning ' fname])
    fclose(fid);
end
c       = fwrite(fid,mn,'int32');
fclose(fid);
status =1;