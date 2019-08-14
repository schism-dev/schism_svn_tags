function dataNew = filter_depth_q(dataOld,depthOld,depthNew)
%dataNew = filter_depth_q(dataOld,depthOld,depthNew)
% interpolates in depth
% quick version, for unique depthOld, and scalar depthNew
%INPUT
% dataOld	- [nTimes, nLevs]
% depthOld	- [nTimes, nLevs] or [1, nLevs]
% depthNew	- [nTimes, nLevs2] or [1, nLevs2]
%
%OUPUT
% dataNew	- [nTimes, nLevs2]
%
% Sergey Frolov
% July 2006, SF, added quick version

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

function dataNew = interpLoop(dataOld,depthOld,depthNew)
  if size(depthNew,2)==1
    dataNew = interpLoopScalar(dataOld,depthOld,depthNew);
  else
    dataNew     = interpLoopVector(dataOld,depthOld,depthNew);
  end


function dataNew = interpLoopScalar(dataOld,depthOld,depthNew)
  %interploate in depth where both depthOld and depthNew are arrays
  %number of rows in depthOld and depthNew == num rows in dataOld
  % this is when depthNew(i,:) is a scalar (single new depth)

  dataNew = nan(size(dataOld,1),size(depthNew,2));
  dataOld=dataOld';
  depthOld=depthOld';

  %loop over datapoints
 for i =1:size(dataOld,2)
    %make sure depthOld(i,:) are unique
    idxNNan=find(~isnan(depthOld(:,i)));
    [junk, idxUnique] = unique(depthOld(idxNNan,i));
    idxUnique=idxNNan(idxUnique);
    if (depthOld(idxUnique(1),i)-depthOld(idxUnique(end),i)); %can't interpolate if first and last layer are hte same
      dataNew(i,1) = interp1qq(depthOld(idxUnique,i),dataOld(idxUnique,i),depthNew(i));
    else
      dataNew(i,1) = nan;
    end
  end

function dataNew     = interpLoopVector(dataOld,depthOld,depthNew)
  %interploate in depth where both depthOld and depthNew are arrays
  %number of rows in depthOld and depthNew == num rows in dataOld
  % this is when depthNew(i,:) is a vector (multiple new depth)

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



function yi=interp1qq(x,y,xi)
%adopted from interp1q for scalar xi
%keyboard
  r = max(find(x <= xi));
  r(xi==x(end)) = length(x)-1;
  if isempty(r) | (r<=0) | (r>=length(x)) 
    yi = NaN(1,size(y,2),superiorfloat(x,y,xi));
  else 
    u = (xi-x(r))./(x(r+1)-x(r)); 
    yi=y(r,:)+(y(r+1,:)-y(r,:)).*u(:,ones(1,size(y,2))); 
  end


