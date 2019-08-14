function fid = fort_addColumnsFid(fname, columns2add,fid)
%fid = fort_add(fname, column2add,[fid])
%adds one more column to the fort file
%makes shure that the size in the file headre is not exceeded
%
%sergey frolov august 2004

if nargin<3
  fid=0;
end
if ~fid
  fid     = fopen(fname,'a+');
end
fseek(fid,0,'bof');
[mn,c]  = fread(fid,2,'int32');     %size m by n decleared in the file header
mnAdd   = size(columns2add);        %size mn of columns to add
fseek(fid,0,'eof');
eofPos  = ftell(fid);
nInFile = (eofPos-4*2)/(4*mn(1));
if (floor(nInFile)-nInFile)~=0
    error(['Cant add a column to the file ' fname ' file has non intger number of columns ' nInFile] )
elseif mn(1)~=mnAdd(1) 
   error(['Cant add a column to the file ' fname ' number of rows in file and new columns doesnt agree'] )
elseif mnAdd(2) + nInFile > mn(2)
   error(['Cant add a column to the file ' fname ' total number of columns ' mnAdd(2) + nInFile ' exceeds n in the file header'] ) 
else
   c   = fwrite(fid,columns2add,'float32');
   status =1;
end

if nargout==0
  fclose(fid);
  fid=0;
end
