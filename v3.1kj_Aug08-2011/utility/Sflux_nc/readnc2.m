%Read in our sflux netcdf files into matlab and plot; modify as appro.
clear all; close all;
%scrsz = get(0,'ScreenSize'); %screen size
%Dimension for arrays reversed from ncdump (FORTRAN convention)

CB_bnd=load('CB_bnd.xy'); %load domain bnd

%NARR files
fill_in=1.e9; %junk value from nc files
avi_out = avifile('sflux.avi');
for i=9:12 %i=9 => Sept. 16
  char=sprintf('%3.3d',i);
  filen=strcat('sflux_air_2.',char,'.nc');
  ncid0 = netcdf.open(filen,'NC_NOWRITE');
  vid=netcdf.inqVarID(ncid0,'time'); %input var./array name
  time_narr= netcdf.getVar(ncid0, vid);
  vid=netcdf.inqVarID(ncid0,'lon');
  lon_narr = netcdf.getVar(ncid0, vid); 
  vid=netcdf.inqVarID(ncid0,'lat');
  lat_narr = netcdf.getVar(ncid0, vid); 

  vid=netcdf.inqVarID(ncid0,'uwind');
  %Time dimension is last index in uwind_narr etc.
  uwind_narr= netcdf.getVar(ncid0, vid); 
  uwind_narr(find(abs(uwind_narr)>fill_in))=nan;
  vid=netcdf.inqVarID(ncid0,'vwind');
  vwind_narr= netcdf.getVar(ncid0, vid); 
  vwind_narr(find(abs(vwind_narr)>fill_in))=nan;
  vid=netcdf.inqVarID(ncid0,'prmsl');
  pres_narr= netcdf.getVar(ncid0, vid); 
  pres_narr(find(abs(pres_narr)>fill_in))=nan;
  vid=netcdf.inqVarID(ncid0,'stmp');
  airt_narr= netcdf.getVar(ncid0, vid); %ait temp.
  airt_narr(find(abs(airt_narr)>fill_in))=nan;
  vid=netcdf.inqVarID(ncid0,'spfh'); %humidity
  spfh_narr= netcdf.getVar(ncid0, vid); %specific humidity
  spfh_narr(find(abs(spfh_narr)>fill_in))=nan;

  disp(strcat('Done reading: ',filen));

  for j=1:length(time_narr)
    %Time stamp
    date_narr=strcat('Sept ',num2str(i+7),': ',num2str(time_narr(j)*24),' hour UTC');
    %Compute bound
    xmin=min(min(lon_narr)); 
    ymin=min(min(lat_narr));
    xmax=max(max(lon_narr));
    ymax=max(max(lat_narr));

    hold on;
    plot(CB_bnd(:,1),CB_bnd(:,2),'k'); %plot bnd

    %Scalar plot
%    contour(lon_narr,lat_narr,pres_narr(:,:,j));
%    caxis([9.7e4 10.2e4]);
%    colorbar;

    %Subsample
    ind1=10:10:size(lon_narr,1);
    ind2=10:10:size(lon_narr,2);

    %Vel. plot
    scale=0.05; %scale to fit
    quiver(lon_narr(ind1,ind2),lat_narr(ind1,ind2),...
    uwind_narr(ind1,ind2,j)*scale,vwind_narr(ind1,ind2,j)*scale,0,'b');
    quiver(-77.9,39.9,5*scale,0.,0,'r');
    text(-77.9,39.8,'5 m/s');
    text(-76.9,39.8,date_narr);
    %Warning: avi does not like changing in xmin etc; use fixed scale
%    axis([xmin-0.1 xmax+0.1 ymin-0.1 ymax+0.1]);
    axis([-81 -72.6 33.32 40.45]);
    xlabel('Lon'); ylabel('Lat');

    frame = getframe(gca);
    avi_out=addframe(avi_out,frame);
    clf; %clear figure
  end %j
end %for all nc files
avi_out=close(avi_out);
