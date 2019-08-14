function h=ob_plot_sb_tr(ob, eta, C,maxmin)
%h=ob_plot_sb_tr(ob, eta, ob_st,[maxmin])
%plots transect for s-coordinate grid
% ob defines mapping for a transect
% elevation at location of transaction nodes
% ob_st - is already extracted transect
% maxmin - optional scaling of the color axis
% h - handle to the pcolor axis
%
%Sergey Frolov june 2005


%compute distance along the trasect
x=zeros(length(ob.xy.x),1);
for i =2:length(ob.xy.x)
    x(i)=x(i-1)+norm(...
        [ob.xy.x(i)-ob.xy.x(i-1),...
        ob.xy.y(i)-ob.xy.y(i-1)]);
end
x=x/1e3;


if min(size(C))==1	%ploting 2d scalar (like elevation)
  h = plot(x,C);
  ylabel('\eta m')
  xlabel('distance km')
  set(gca,'FontSize',14);
%  axis([min(x) max(x) -30 2]);
  if nargin ==4
      axis([min(x) max(x) maxmin]);
  end
else			%ploting 3d trasect (like salt, temp, ets)
  dp	= ob.xy.H*ob.gr.hgrid.depth;
  Y = sb_computeZlevels(dp, eta, ob.gr.vgrid)';
  [n m] = size(C);
  X = repmat(x',n,1);
  faces = nan*ones(n-1*m-1,4);
  k = 1;
  for i = 1:m-1
    for j = 1:n-1
      faces(k,:) = [(i-1)*n+j, (i-1)*n+j+1,  i*n+j+1, i*n+j,];
      k=k+1;
    end
  end

  %get curent axes, if no ases exist create them, otherwise make it comaptible with subplot
  ha=gca;
  xl=get(ha,'xLim');      %get curent limits
  yl=get(ha,'yLim');
  cla                     %clear axes so we can use the patch function

  patch('Faces',faces,'Vertices',[X(:) Y(:)], 'FaceVertexCData', C(:),'FaceColor','interp','EdgeColor','none')

  if (xl(1)==0|xl(2)==1|yl(1)==0|yl(2)==1)
    %axis equal
  else
    axis([xl yl])
  end

  xlabel('distance km')
  ylabel('depth m')
%  set(gca,'FontSize',14);
%  axis([min(x) max(x) -30 2]);
  if nargin ==4
      caxis(maxmin)
  end
  h2=colorbar;
  set(h2,'FontSize',14);
end

