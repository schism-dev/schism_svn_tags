function sb_addTimeStep(fname, t, it, eta, data, fid)
% sb_addTimeStep(fname, t,it,eta,data, [fid])
% add time steps to selfe  binary file fname
% fname - file name
% t 	- time (float)
% it 	- iteration number
% eta	- elevation vector (1:np,1)
% data 	- data from the file [1:np*nLevels,1:numVectComp]
% fid	- optionaly you can pass a valid fid to write to
%
% SF, november 2005,

%tic

if nargin < 6
 fid = fopen(fname,'a');
end

c = fwrite(fid,t,'float32');
c = fwrite(fid,it,'int32');
c = fwrite(fid, eta, 'float32');
c = fwrite(fid, data', 'float32')';

if nargin < 6
  fid=fclose(fid);
end


