%Same as readnc2 but also output time series at a given point (nearest neighbor interp)
clear all; close all;
%scrsz = get(0,'ScreenSize'); %screen size
%Dimension for arrays reversed from ncdump (FORTRAN convention)

%Input point loc
xin=-9.55956300e+001;
yin=2.91482159e+001;

CB_bnd=load('CB_bnd.xy'); %load domain bnd

%NARR files
fill_in=1.e9; %junk value from nc files
avi_out = avifile('sflux.avi');
iout=0; %counter for output
for i=1:12 %stack # for nc files
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

  %Find nearest node for interpolation
  dist2=(lon_narr-xin).^2+(lat_narr-yin).^2;
  [tmp,irow0]=min(dist2);
  [tmp2,icol]=min(tmp);
  irow=irow0(icol); 
%  disp('Nearest pt='); [lon_narr(irow,icol) lat_narr(irow,icol)]

  for j=1:length(time_narr)
    %Time stamp for plot - modify as appropr.
    date_narr=strcat('Nov ',num2str(i+15),': ',num2str(time_narr(j)*24),' hour UTC');
    iout=iout+1;
    time_out(iout)=(i-1)+time_narr(j); %in days (UTC) from _001.nc
    uwind_out(iout)=uwind_narr(irow,icol,j);
    vwind_out(iout)=vwind_narr(irow,icol,j);

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
    ind1=1:1:size(lon_narr,1);
    ind2=1:1:size(lon_narr,2);

    %Vel. plot
    scale=0.02; %scale to fit
    quiver(lon_narr(ind1,ind2),lat_narr(ind1,ind2),...
    uwind_narr(ind1,ind2,j)*scale,vwind_narr(ind1,ind2,j)*scale,0,'b');
    quiver(-95.35,29.112,5*scale,0.,0,'r');
    text(-95.35,29.13,'5 m/s');
    text(-95.35,29.,date_narr);
    %Warning: avi does not like changes in xmin etc; use fixed scale
%    axis([xmin-0.1 xmax+0.1 ymin-0.1 ymax+0.1]);
    axis([-95.63 -95.17 28.748 29.15]);
    xlabel('Lon'); ylabel('Lat');

    frame = getframe(gca);
    avi_out=addframe(avi_out,frame);
    clf; %clear figure
  end %j

  clear time_narr lon_narr lat_narr uwind_narr vwind_narr pres_narr airt_narr spfh_narr dist2 irow icol tmp tmp2;
end %for all nc files
avi_out=close(avi_out);

%Plot time series
narr=load('wind.uv_from_Nov16_2008'); %time,u,v
figure(2);
plot(time_out,uwind_out,'r',time_out,vwind_out,'k',narr(:,1),narr(:,2),'b',...
     narr(:,1),narr(:,3),'g');
legend('Uwind (m/s)','Vwind');
xlabel('Time (days UTC)');
