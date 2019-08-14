%Read in our sflux netcdf files into matlab and plot; modify as appro.
% e.g., the start and end stack #; time stamp for plot; sub-sampling freq.
% in plot
clear all; close all;
%scrsz = get(0,'ScreenSize'); %screen size
%Dimension for arrays reversed from ncdump (FORTRAN convention)

CB_bnd=load('CB_bnd.xy'); %load domain bnd

%NARR files
fill_in=1.e9; %junk value from nc files
avi_out = avifile('sflux.avi');
for i=1:20 %stack # for nc files
  char=sprintf('%3.3d',i);
  filen=strcat('sflux_air_1.',char,'.nc');
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
    %Time stamp for plot - modify as appropr.
    date_narr=strcat('Sept ',num2str(i),': ',num2str(time_narr(j)*24),' hour UTC');
    %Compute bound
    xmin=min(min(lon_narr)); 
    ymin=min(min(lat_narr));
    xmax=max(max(lon_narr));
    ymax=max(max(lat_narr));

    %plot domain bnd
    hold on;
    plot(CB_bnd(:,1),CB_bnd(:,2),'k'); 

    %Scalar plot
%    contour(lon_narr,lat_narr,pres_narr(:,:,j));
%    caxis([9.7e4 10.2e4]);
%    colorbar;

    %Subsample
    ind1=1:1:size(lon_narr,1);
    ind2=1:1:size(lon_narr,2);

    %Vel. plot
    scale=0.2; %scale to fit
    quiver(lon_narr(ind1,ind2),lat_narr(ind1,ind2),...
    uwind_narr(ind1,ind2,j)*scale,vwind_narr(ind1,ind2,j)*scale,0,'b');
    quiver(-74.7,5,5*scale,0.,0,'r');
    text(-74.7,5.2,'5 m/s');

    %Atmos. pressure
    pcolor(lon_narr(ind1,ind2),lat_narr(ind1,ind2),double(pres_narr(ind1,ind2,j)/1.e5));
    colorbar;

    text(-56,8,date_narr);
    %Warning: avi does not like changes in xmin etc; use fixed scale
%    axis([xmin-0.1 xmax+0.1 ymin-0.1 ymax+0.1]);
    %axis([-123.01 -121.4 37.38 38.25]);
    axis([-74.66 -53.28 3.73 24.5]);
    xlabel('Lon'); ylabel('Lat');
    title('Atmos. pressure scaled by 1.e5');

    frame = getframe(gca);
    avi_out=addframe(avi_out,frame);
    clf; %clear figure
  end %j

  clear time_narr lon_narr lat_narr uwind_narr vwind_narr pres_narr airt_narr spfh_narr;
end %for all nc files
avi_out=close(avi_out);
