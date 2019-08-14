function ebHeader   = eb_changeFname(ebHeader,fname)
%ebHeader   = eb_changeFname(ebHeader,fname)
%upadtes file name (fname) in elcirc header ebHeader
%only limited checking of the errors
%
%Sergey Frolov,     August 2004

if exist(fname,'file')~=2
    error(['file' fname 'doesnt exist']);
end

ebHeader.fname = fname;
