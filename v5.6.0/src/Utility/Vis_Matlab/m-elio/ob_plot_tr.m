function h=ob_plot_tr(ob_tr, ob_st,maxmin)
%h=ob_plot_tr(ob_tr, ob_st,[maxmin])
%plots transect
% ob_tr defines mapping for a transect
% ob_st - is already extracted transect
% maxmin - optional scaling of the color axis
% h - handle to the pcolor axis
%
%Sergey Frolov dec 2004



x=zeros(length(ob_tr.xy.x),1);
for i =2:length(ob_tr.xy.x)
    x(i)=x(i-1)+norm(...
        [ob_tr.xy.x(i)-ob_tr.xy.x(i-1),...
        ob_tr.xy.y(i)-ob_tr.xy.y(i-1)]);
end
x=x/1e3;
y=ob_tr.gr.vgrid.zLevel-ob_tr.gr.vgrid.zMsl;

if min(size(ob_st))==1	%ploting 2d scalar (like elevation)
  h = plot(x,ob_st);
  ylabel('\eta m')
  xlabel('distance km')
  set(gca,'FontSize',14);
%  axis([min(x) max(x) -30 2]);
  if nargin ==3
      axis([min(x) max(x) maxmin]);
  end
else			%ploting 3d trasect (like salt, temp, ets)
  h=pcolor(x,y,ob_st');
  xlabel('distance km')
  ylabel('depth m')
  set(gca,'FontSize',14);
  set(h,'LineStyle','none','FaceColor','interp');
  axis([min(x) max(x) -30 2]);
  if nargin ==3
      caxis(maxmin)
  end
  h2=colorbar;
  set(h2,'FontSize',14);
end
