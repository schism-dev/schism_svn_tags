function [data,ts] = eb_readTimeStep2(h,n,fn)
% [data,ts] = eb_readTimeStep2(h,n)
% read time steps n=[] from elcirc binary file h.fname
% h 	- header structure as returned by eb_readHeader
% n 	- array of timestep numbers 
% data 	- data from the file [dataVector,vectorComponent_k,time_iteration]
%	in before matlab v14 output is in double precision othervise in single
% ts 	- some info on time step
% sergey frolov march 8, 2004
% SF, may 2005, migrated computeStepSize computeStepIdx to eb_readHeader
%	also added switch for using new functionality in matlab v 14

%tic

if n>h.nSteps
    error('time step out of range')
end
if ~exist('fn','var')
  fn = h.fname;
end

%if ~h.fid	%if h is not modified in this function then it is passed by reference => better speed
    fid   = fopen(fn);
%else 
%    fid   = h.fid;
%end

%read data, using diferent functions depending on the version of matlab
if str2num(version('-release')) < 14
  data=zeros([h.gridSize h.flagSv length(n)]);
  for i=1:length(n)
  	[t,d]	=readTs_v13(h,n(i),fid);
        data(:,:,i)=d;
        ts{i}	=t;
  end
else
  data=zeros([h.gridSize h.flagSv length(n)],'single');
  for i=1:length(n)
        [t,d]   =readTs_v14(h,n(i),fid);
        data(:,:,i)=d;
        ts{i}   =t;
  end
end

%close file
if ~h.fid
    fid=fclose(fid);
end

% disp(['Read ' num2str(length(n)) ' time steps in ' num2str(toc)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%private functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ts,data]=readTs_v13(h,n,fid)
%read time step return results in double precision
    fseek(fid, h.dataStartPos + (h.stepSize)*(n-1), -1);
    
    ts.t        = fread(fid,1,'float32');
    ts.tit      = fread(fid,1,'int32');
    ts.sidx     = uint8(fread(fid, h.hgrid.np, 'int32')); 
    data    	= fread(fid, [h.flagSv h.gridSize], 'float32')';            

function [ts,data]=readTs_v14(h,n,fid)
%read time step return results in single precision
%should use memory maps in future
    fseek(fid, h.dataStartPos + (h.stepSize)*(n-1), -1);

    ts.t        = fread(fid,1,'float32=>float32');
    ts.tit      = fread(fid,1,'int32=>int32');
    ts.sidx     = uint8(fread(fid, h.hgrid.np, 'int32=>int32'));
    data        = fread(fid, [h.flagSv h.gridSize], 'float32=>float32')';



