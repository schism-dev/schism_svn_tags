function [ startTime ] = eb_getStartTime( fname )
% [ startTime ] = eb_getStartTime( fname )
% gets start time from the header of the elcirc binary file
%
% Sergey Frolov mar 08, 2005

fid          = fopen(fname);

dataFormat   = char(fread(fid,48,'uchar')');
version      = char(fread(fid,48,'uchar')');
startTime    = char(fread(fid,48,'uchar')');

fclose(fid);
    
