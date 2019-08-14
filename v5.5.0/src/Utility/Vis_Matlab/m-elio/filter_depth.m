function dataNew = filter_depth(dataOld,depthOld,depthNew)
%dataNew = filter_depth(dataOld,depthOld,depthNew)
% interpolates in depth
%INPUT
% dataOld	- [nTimes, nLevs]
% depthOld	- [nTimes, nLevs] or [1, nLevs]
% depthNew	- [nTimes, nLevs2] or [1, nLevs2]
%
%OUPUT
% dataNew	- [nTimes, nLevs2]
%
% Sergey Frolov

if size(depthOld,2)~=size(dataOld,2)
  error('size(depthOld,2) should be == to size(dataOld,2)')
end

if (size(depthOld,1)==1)&(size(depthNew,1)==1)
  dataNew	= interp1(depthOld,dataOld,depthNew);
elseif (size(depthOld,1)~=1)&(size(depthNew,1)==1)
  if size(depthOld,1)~=size(dataOld,1)
    error('size(depthOld,1) should be == to size(dataOld,1)')
  else
    depthNew    = repmat(depthNew,size(dataOld,1),1);
    dataNew	= interpLoop(dataOld,depthOld,depthNew);
  end
elseif (size(depthOld,1)==1)&(size(depthNew,1)~=1)
  if size(depthNew,1)~=size(dataOld,1)
    error('size(depthNew,1) should be == to size(dataOld,1)')
  else
    depthOld    = repmat(depthOld,size(dataOld,1),1);
    dataNew     = interpLoop(dataOld,depthOld,depthNew);
  end
elseif (size(depthOld,1)~=1)&(size(depthNew,1)~=1)
  if (size(depthOld,1)~=size(dataOld,1))|size(depthNew,1)~=size(dataOld,1)
    error('size(depthNew,1) and size(depthOld,1) should be == to size(dataOld,1)')
  else
%    depthNew    = repmat(depthNew,size(dataOld,1),1);
    dataNew     = interpLoop(dataOld,depthOld,depthNew);
  end
else
  error('confused about size of your input arrays')
end


function dataNew     = interpLoop(dataOld,depthOld,depthNew)
%interploate in depth where both depthOld and depthNew are arrays
%number of rows in depthOld and depthNew == num rows in dataOld

dataNew = nan*zeros(size(dataOld,1),size(depthNew,2));
%loop over datapoints
warning('off','MATLAB:interp1:NaNinY')
for i =1:size(dataOld,1)
%  i = idx(j);
	%make sure depthOld(i,:) are unique
  idxNNan=find(~isnan(depthOld(i,:)));
  [junk, idxUnique] = unique(depthOld(i,idxNNan));
  idxUnique=idxNNan(idxUnique);
  if length(idxUnique)>=2
    dataNew(i,:) = interp1(depthOld(i,idxUnique),dataOld(i,idxUnique),depthNew(i,:));
  else
    dataNew(i,:) = nan*depthNew(i,:);
  end
end
warning('on','MATLAB:interp1:NaNinY')



