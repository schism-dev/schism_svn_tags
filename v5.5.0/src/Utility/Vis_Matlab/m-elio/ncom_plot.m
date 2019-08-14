function h=ncom_plot(hgr,var,lims)

h=pcolor(hgr.lon,hgr.lat,var);
set(h,'edgeColor','none')

if exist('lims','var')
  clims(lims)
end
colorbar

