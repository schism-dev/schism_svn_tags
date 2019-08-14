function status = fort_add(fname, column2add)
%status = fort_add(fname, column2add)
%adds one more column to the fort file
%makes shure that the size in the file headre is not exceeded
%
%sergey frolov august 2004

fid     = fopen(fname,'a+');
[mn,c]  = fread(fid,2,'int32');  %size m by n
fseek(fid,0,'eof');
eofPos  = ftell(fid);

if eofPos>=4*(2+mn(1)*mn(2))
    error(['Cant add a column to the file ' fname ' column number exceeds n in the file header'] )
    status=0;
elseif size(column2add,1)~=mn(1)
    error(['Cant add a column to the file ' fname ' column length ~= m in the file header'] )
    status=0;    
else
    c   = fwrite(fid,column2add,'float32');
    status =1;
end
fclose(fid);