function [pobj]=gr_plot(gr,c,maxmin)
%function [h]=gr_plot(gr,c,[maxmin])
% displays individula grid of gr=readGrid structure
% gr is one of the gr.hgrid, gr.ecenters, ....
% c - index color matrix for each node
% maxmin[max min] -maximum and minimum of the color stretch
% Example: gr_plot(gr.ecenters,gr.ecenters.depth)
%
% Sergey Frolov,   March 2004

r=max(c)-min(c);
ncol=200;

if nargin ==3
   c( find(c < maxmin(1)) )= maxmin(1);
   c( find(c > maxmin(2)) )= maxmin(2);   
end

%get curent axes, if no ases exist create them, otherwise make it comaptible with subplot
ha=gca;
xl=get(ha,'xLim');	%get curent limits
yl=get(ha,'yLim');
cla			%clear axes so we can use the patch function

% do tris
pobj(1) = patch('Faces',gr.elem(find(gr.elem(:,2)==3),3:end-1),'Vertices',gr.nodes(:,2:end-1),'FaceVertexCData',c,'FaceColor','interp','EdgeColor','none');
% do quads
pobj(2) = patch('Faces',gr.elem(find(gr.elem(:,2)==4),3:end),'Vertices',gr.nodes(:,2:end-1),'FaceVertexCData',c,'FaceColor','interp','EdgeColor','none');
%set(gca,'FontSize',14);

if isfield(gr,'bndLine')
  if ~isempty(gr.bndLine.X)
    line(gr.bndLine.X,gr.bndLine.Y,'color','k')
  end
end

% if xl=yl=[0 1] this axes are new and we don;t want to scale them
% otherwise use existing xl yl so you don't need to rescale graph evry time you plot
if (xl(1)==0|xl(2)==1|yl(1)==0|yl(2)==1)
  axis equal
else
  axis([xl yl])
end

%if limits are provided use user specified color axes
if nargin ==3
  caxis(maxmin)
else
  caxis auto
end

%colormap(jet(ncol))
%h=colorbar;
%set(h,'FontSize',14);
%axis equal
