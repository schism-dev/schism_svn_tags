function hgrid=gr_parseBnd(hgrid)

if ~isfield(hgrid,'eofLines')
  error('No boundary information')
end
if isempty(hgrid.eofLines)
  error('No boundary information')
end

X=[];
Y=[];

%number of open boundaries
idxStart=1;
nb=sscanf(hgrid.eofLines{1},'%d',1);
idxStart=idxStart+1;
for i=1:nb
  nnodes=sscanf(hgrid.eofLines{idxStart+1},'%d',1);
  idxStart=idxStart+1;
  nodes=nan(nnodes,1);
  for j=1:nnodes
    nodes(j)=sscanf(hgrid.eofLines{idxStart+j},'%d',1);
  end
  idxStart=idxStart+nnodes;

  x=hgrid.x(nodes);
  y=hgrid.y(nodes);

  X_=[x(1:end-1)';x(2:end)'];
  Y_=[y(1:end-1)';y(2:end)'];

  X=[X X_];
  Y=[Y Y_];

end
idxStart=idxStart+1;

%number of land boundaries
nb=sscanf(hgrid.eofLines{idxStart},'%d',1);
idxStart=idxStart+1;
for i=1:nb
  tmp=sscanf(hgrid.eofLines{idxStart+1},'%d',2);
  nnodes=tmp(1);islandFlag=tmp(2);
  idxStart=idxStart+1;
  nodes=nan(nnodes,1);
  for j=1:nnodes
    nodes(j)=sscanf(hgrid.eofLines{idxStart+j},'%d',1);
  end
  idxStart=idxStart+nnodes;

  x=hgrid.x(nodes);
  y=hgrid.y(nodes);

  if islandFlag
    X_=[x([1:end-1 1])';x([2:end end])'];
    Y_=[y([1:end-1 1])';y([2:end end])'];
  else %land boundary
    X_=[x(1:end-1)';x(2:end)'];
    Y_=[y(1:end-1)';y(2:end)'];
  end

  X=[X X_];
  Y=[Y Y_];

end

hgrid.bndLine.X=X;
hgrid.bndLine.Y=Y;

%line(X,Y,'Color','k')


