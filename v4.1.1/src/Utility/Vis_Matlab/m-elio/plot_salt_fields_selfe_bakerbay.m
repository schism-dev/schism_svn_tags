% plot_salt_fields_selfe.m  11/19/2008 Parker MacCready
%
% This plots a figure of model surface salinity, a salinity section, and an
% observed salinity section
%
% The intention is to hand this to Antonio's group and they will add in
% panels from their model (SELFE)
%

clear
%close all
path(path,'./m-elio');
%--Coastline Information
load '/home/workspace/ccalmr45/choj/for_grant/barb_sal_maps/c28142.dat';
c1=c28142;
x1=c1(:,1);y1=c1(:,2);

base_out_dir ='/home/workspace/ccalmr45/choj/selfe/2009-27-37-301a/';
%base_out_dir ='/home/workspace/ccalmr45/choj/selfe/for_scott/bakerbay_refined/';
start_day = datenum(2009, 1, 1, 0, 0, 0);

%-- This section will be changed for each day and each forecast.

curr_run_links_dir = '/home/workspace/ccalmr45/choj/selfe/2009-27-37-301b/2009-32-301b/run';
curr_run_links_dir2 = '/home/workspace/ccalmr45/choj/selfe/2009-27-37-301c/2009-32-301c/run';
%curr_run_links_dir = '/home/workspace/ccalmr45/choj/selfe/2009-10-26-300a';
%curr_run_links_dir = '/home/workspace/ccalmr45/choj/selfe/for_scott/bakerbay_refined';
yr=2009;imn=8;idy=6; % variables: today's date
snday=1;enday=7;ifg=1;
ist=4;intv=4;ied=96; % ist=start;intv=interval; ied=end; time loop: it=1:1:96--> every 15min; it=4:4:96--> hourly
if (ifg==1) 
    aviobj=avifile('ani_sal4bakerbay_sbs_comp2_32w.avi','fps',1);
end
vbot=18; % vbot=18 for fca200nb and dev;
pt1c=12920;pt2c=12967;pt3c=15016;pt4c=15702;
pt1=13326;pt2=14853;pt3=22434;pt4=21807;
%-------------------------------------------------------------

curr_day=datenum(yr, imn, idy, 0, 0, 0);

eval(['!ln -sf ' curr_run_links_dir '/hgrid.gr3 hgrid2.gr3'])
eval(['!ln -sf ' curr_run_links_dir '/hgrid.ll hgrid2.ll'])
eval(['!ln -sf ' curr_run_links_dir2 '/hgrid.gr3 hgrid2c.gr3'])
eval(['!ln -sf ' curr_run_links_dir2 '/hgrid.ll hgrid2c.ll'])
% prepare fields - salt & veolocity
gr.hgrid=gr_readHGrid('hgrid2.gr3');
gr2.hgrid=gr_readHGrid('hgrid2.ll');
fg = hgrid2fg('hgrid2.ll');
grc.hgrid=gr_readHGrid('hgrid2c.gr3');
gr2c.hgrid=gr_readHGrid('hgrid2c.ll');
fgc = hgrid2fg('hgrid2c.ll');
% set some default values
aa = [-124.06 -123.92 46.23 46.33]; % map region
dar = [1/cos(pi*46/180) 1 1]; % Cartesian scaling
clevs = [-200 -200]; % bathymetry contours
slo = 0; shi = 34; % salinity color range
sclevs = [30:.5:34]; % salinity contour levels
sclevs2 = [10:2:20]; % salinity contour levels
fs1 = 14; % fontsize
dlat0 = 0.001; % spacing for regridded velocity (degrees of latitude)
[LON0,LAT0] = meshgrid([aa(1):dlat0*dar(1):aa(2)],[aa(3):dlat0:aa(4)]);
z0 = griddata(fg.x,fg.y,-gr.hgrid.depth,LON0,LAT0,'cubic');
z0c = griddata(fgc.x,fgc.y,-grc.hgrid.depth,LON0,LAT0,'cubic');
dlat2 = 0.05; % spacing for regridded velocity (degrees of latitude)
[LON2,LAT2] = meshgrid([aa(1):dlat2*dar(1):aa(2)],[aa(3):dlat2:aa(4)]);
dlat = 0.1; % spacing for regridded velocity (degrees of latitude)
[LON,LAT] = meshgrid([aa(1):dlat*dar(1):aa(2)],[aa(3):dlat:aa(4)]);

nnn=0;
for j=snday:enday %j
 eval(['!ln -sf ' curr_run_links_dir '/' sprintf('%d',j) '_hvel.64 ' sprintf('%d',j) '_hvel.64'])
 eval(['!ln -sf ' curr_run_links_dir '/' sprintf('%d',j) '_salt.63 ' sprintf('%d',j) '_salt.63'])
 eval(['!ln -sf ' curr_run_links_dir '/' sprintf('%d',j) '_elev.61 ' sprintf('%d',j) '_elev.61'])
 h=sz_readHeader([sprintf('%d_',j) 'salt.63']);
 gr.vgrid=h.vgrid;
 hv=sz_readHeader([sprintf('%d_',j) 'hvel.64']);
 
 eval(['!ln -sf ' curr_run_links_dir2 '/' sprintf('%d',j) '_hvel.64 ' sprintf('%d',j) '_hvelc.64'])
 eval(['!ln -sf ' curr_run_links_dir2 '/' sprintf('%d',j) '_salt.63 ' sprintf('%d',j) '_saltc.63'])
 eval(['!ln -sf ' curr_run_links_dir2 '/' sprintf('%d',j) '_elev.61 ' sprintf('%d',j) '_elevc.61'])
 hc=sz_readHeader([sprintf('%d_',j) 'saltc.63']);
 grc.vgrid=hc.vgrid;
 hvc=sz_readHeader([sprintf('%d_',j) 'hvelc.64']);
 idy=idy
 if (idy > 31)
   idy=idy-31;
   imn=imn+1;
 end
 for it=ist:intv:ied  % time loop: it=1:96--> every 15min; it=4:4:96--> hourly 
    nnn=nnn+1;
    if (ifg==1)
     figure(100);
     set(gcf,'position',[100 300 1100 600],'Color','w');
    end
    ihr=floor(it/4);
    imt=(it/4-ihr)*60;

    if (ihr==24)
        ihr=0;
        idy=idy+1;
    end
    if (idy > 31)
        idy=idy-31;
        imn=imn+1;
    end
    if ifg==1
    ctime=sprintf('%d-%02d-%d %02d:%02d %s', yr, imn, idy, ihr, imt, 'PST');

    set(0,'defaultaxesfontsize',fs1);
    set(0,'defaulttextfontsize',fs1);
    set(0,'defaultaxesfontname','Times');
    set(0,'defaulttextfontname','Times');
    colormap(jet(14)); % make a steppy colormap
    end
    % SELFE Surface salinity map with velocity vectors
    %
    [d ts]=sz_readTimeStep(h,it); %step 64 is 1600 hrs PST (= 0000 GMT)
    [u] = map_sz2hts(h,d,1);
    [s0] = u(:,end); % End is the surface layer
    [s1] = u(:,vbot); % vbot is the bottom layer
    s0(s0<0) = NaN;
    s1(s1<0) = NaN;
    sal=s0;
    salb=s1;

    [dc ts]=sz_readTimeStep(hc,it); %step 64 is 1600 hrs PST (= 0000 GMT)
    [uc] = map_sz2hts(hc,dc,1);
    [s0c] = uc(:,end); % End is the surface layer
    [s1c] = uc(:,vbot); % vbot is the bottom layer
    s0c(s0c<0) = NaN;
    s1c(s1c<0) = NaN;
    salc=s0c;
    salbc=s1c;
%    s0 = griddata(fg.x,fg.y,s0,LON0,LAT0,'cubic');
%    s1 = griddata(fg.x,fg.y,s1,LON0,LAT0,'cubic');
%    s0(z0>-5) = NaN;
    % prepare fields - velocity
    if 0>1
    [d ts]=sz_readTimeStep(hv,it); %step 64 is 1600 hrs PST (= 0000 GMT)
    u = map_sz2hts(hv,d(:,1),1);
    v = map_sz2hts(hv,d(:,2),1);
    u = u(:,end);
    v = v(:,end);
    utop = griddata(fg.x,fg.y,u,LON0,LAT0,'cubic');
    vtop = griddata(fg.x,fg.y,v,LON0,LAT0,'cubic');
    utop0=utop;vtop0=vtop;
    utop0(z0>-20)=NaN;
    vtop0(z0>-20)=NaN;
    utop(z0>-5)=NaN;
    vtop(z0>-5)=NaN;
    uu = interp2(LON0,LAT0,utop0,LON,LAT);
    vv = interp2(LON0,LAT0,vtop0,LON,LAT);
    uu2 = interp2(LON0,LAT0,utop,LON2,LAT2);
    vv2 = interp2(LON0,LAT0,vtop,LON2,LAT2);
    end
    if ifg==1
    % and plot them
    subplot(121)
    set(gca,'Position',[0.08 0.08 0.4 0.87])
    gr_plot2(gr2c.hgrid,salbc,[slo shi]);
    caxis([slo shi]);
    colorbar('north')
    axis(aa); set(gca,'dataaspectratio',dar);
    hold on;
    % add lat/lon grid lines
%    for ii=1:10
%    plot([-126.00+0.25*ii -126.00+0.25*ii],[40 50],'w--','LineWidth',0.1);hold on; 
%    end
%    for ii=1:15
%    plot([-130 -120],[45.0+0.25*ii 45.0+0.25*ii],'w--','LineWidth',0.1);hold on; 
%    end
    % add velocity vectors
%    ufact = .1;
%    quiver(LON,LAT,ufact*uu*dar(1),ufact*vv,0,'k');
    % add a velocity scale
%    uscale = 0.5;
%    hold on; quiver(-123.75,45.65,ufact*uscale*dar(1),ufact*0,0,'k');
%    hold on; quiver(-123.75,45.65,ufact*0*dar(1),ufact*uscale,0,'k');
%    text(-123.9,45.6,[num2str(uscale),' m s^{-1}'],'color','k');
    % add salinity contours
%    contour(LON0,LAT0,s0,sclevs,'-k');
    % add bathymetry contours
%    [cc,hh] = contour(LON0,LAT0,z0,[-50 -50],'-w');
%    for ii = 1:length(hh); set(hh(ii),'linewidth',2); end;
%    [cc,hh] = contour(LON0,LAT0,z0,clevs,'-w');
%    for ii = 1:length(hh); set(hh(ii),'linewidth',2); end;
%    [cc,hh] = contour(LON0,LAT0,z0,[-1000 -1000],'-w');
%    for ii = 1:length(hh); set(hh(ii),'linewidth',2); end;
    % add notes and labels
%    title('(a) Bottom Salinity');
    xlabel('Longitude (deg)');
    ylabel('Latitude (deg)');
    % add section line
%    plot([slonw(1) slone(1)],[slat(1) slat(1)],'-k','linewidth',2)
%    plot([slonw(2) slone(2)],[slat(2) slat(2)],'-k','linewidth',2)
%    plot([slonw(3) slone(3)],[slat(3) slat(3)],'-k','linewidth',2)
    % and section text
%    text(-123.87,slat(nn),sect)
    hold on;plot(x1,y1,'k-','LineWidth',2);hold on;grid on;
%    text(-125.48,45.1,'Isobaths = -1000  -200  -50 (m)','FontSize',14,'FontWeight','b','Color','k');hold on;
    text(-123.98,46.312,ctime,'FontSize',14,'FontWeight','b','Color','k');hold on;
%    plot(x2,y2,'k+','LineWidth',3,'MarkerSize',14);hold on;
%%    plot(xx2(nnn),yy2(nnn),'kO','LineWidth',3,'MarkerSize',14);hold on;

    subplot(122)
    set(gca,'Position',[0.55 0.08 0.4 0.87])
%    pcolor(LON0,LAT0,s0);
%    shading interp
    gr_plot2(gr2.hgrid,salb,[slo shi]);
    caxis([slo shi]);
    colorbar('north')
    axis(aa); set(gca,'dataaspectratio',dar);
    hold on;
    % add lat/lon grid lines
%    for ii=1:5
%    plot([-124.5+0.1*ii -124.5+0.1*ii],[40 50],'w--','LineWidth',0.1);hold on; 
%    end
%    for ii=1:7
%    plot([-130 -120],[45.8+0.1*ii 45.8+0.1*ii],'w--','LineWidth',0.1);hold on; 
%    end
    % add velocity vectors
%    ufact = .1;
%    quiver(LON2,LAT2,ufact*uu2*dar(1),ufact*vv2,0,'k');
    % add a velocity scale
%    uscale = 0.5;
%    hold on; quiver(-124.45,46.05,ufact*uscale*dar(1),ufact*0,0,'k');
%    hold on; quiver(-124.45,46.05,ufact*0*dar(1),ufact*uscale,0,'k');
%    text(-124.48,46.00,[num2str(uscale),' m s^{-1}'],'color','k');
    % add salinity contours
%    contour(LON0,LAT0,s0,sclevs,'-k');
    % add bathymetry contours
%    [cc,hh] = contour(LON0,LAT0,z0,[-50 -50],'-w');
%    for ii = 1:length(hh); set(hh(ii),'linewidth',2); end;
%    [cc,hh] = contour(LON0,LAT0,z0,clevs,'-w');
%    for ii = 1:length(hh); set(hh(ii),'linewidth',2); end;
%    [cc,hh] = contour(LON0,LAT0,z0,[-100 -100],'-w');
%    for ii = 1:length(hh); set(hh(ii),'linewidth',2); end;
    % add notes and labels
%    title('(b) Bottom Salinity');
    xlabel('Longitude (deg)');
    ylabel('Latitude (deg)');
    % add section line
%    plot([slonw(1) slone(1)],[slat(1) slat(1)],'-k','linewidth',2)
%    plot([slonw(2) slone(2)],[slat(2) slat(2)],'-k','linewidth',2)
%    plot([slonw(3) slone(3)],[slat(3) slat(3)],'-k','linewidth',2)
    % and section text
%    text(-123.87,slat(nn),sect)
    hold on;plot(x1,y1,'k-','LineWidth',2);hold on;grid on;
    text(-123.98,46.312,ctime,'FontSize',14,'FontWeight','b','Color','k');hold on;
%    plot(xx2(nnn),yy2(nnn),'kO','LineWidth',3,'MarkerSize',14);hold on;
    
    frame=getframe(gcf);
    aviobj=addframe(aviobj,frame);
    end % ifg==1
    %--Get time series of specific locations
    ttm(nnn)=datenum(yr,imn,idy,ihr,imt,0);
    ss1c(nnn)=salc(pt1c);sb1c(nnn)=salbc(pt1c);
    ss2c(nnn)=salc(pt2c);sb2c(nnn)=salbc(pt2c);
    ss3c(nnn)=salc(pt3c);sb3c(nnn)=salbc(pt3c);
    ss4c(nnn)=salc(pt4c);sb4c(nnn)=salbc(pt4c);
    ss1(nnn)=sal(pt1);sb1(nnn)=salb(pt1);
    ss2(nnn)=sal(pt2);sb2(nnn)=salb(pt2);
    ss3(nnn)=sal(pt3);sb3(nnn)=salb(pt3);
    ss4(nnn)=sal(pt4);sb4(nnn)=salb(pt4);
 end % loop: it
end % loop:j
if ifg==1
    aviobj=close(aviobj);
    close all;
end
saldata(:,1)=ttm(:)';
saldata(:,2)=ss1c(:)';
saldata(:,3)=sb1c(:)';
saldata(:,4)=ss2c(:)';
saldata(:,5)=sb2c(:)';
saldata(:,6)=ss3c(:)';
saldata(:,7)=sb3c(:)';
saldata(:,8)=ss4c(:)';
saldata(:,9)=sb4c(:)';
saldata(:,10)=ss1(:)';
saldata(:,11)=sb1(:)';
saldata(:,12)=ss2(:)';
saldata(:,13)=sb2(:)';
saldata(:,14)=ss3(:)';
saldata(:,15)=sb3(:)';
saldata(:,16)=ss4(:)';
saldata(:,17)=sb4(:)';

figure(101)
set(gcf,'position',[100 200 550 700],'Color','w')
subplot(4,1,1)
plot(ttm,ss1c,'b-',ttm,ss1,'r-')
datetick('x','mm/dd','keeplimits','keepticks');
ylabel('S (psu)','FontSize',14);
subplot(4,1,2)
plot(ttm,ss2c,'b-',ttm,ss2,'r-')
datetick('x','mm/dd','keeplimits','keepticks');
ylabel('S (psu)','FontSize',14);
subplot(4,1,3)
plot(ttm,ss3c,'b-',ttm,ss3,'r-')
datetick('x','mm/dd','keeplimits','keepticks');
ylabel('S (psu)','FontSize',14);
subplot(4,1,4)
plot(ttm,ss4c,'b-',ttm,ss4,'r-')
datetick('x','mm/dd','keeplimits','keepticks');
ylabel('S (psu)','FontSize',14);
xlabel('Date in 2009 (PST)','FontSize',14);

figure(102)
set(gcf,'position',[100 200 550 700],'Color','w')
subplot(4,1,1)
plot(ttm,sb1c,'b-',ttm,sb1,'r-')
datetick('x','mm/dd','keeplimits','keepticks');
ylabel('S (psu)','FontSize',14);
subplot(4,1,2)
plot(ttm,sb2c,'b-',ttm,sb2,'r-')
datetick('x','mm/dd','keeplimits','keepticks');
ylabel('S (psu)','FontSize',14);
subplot(4,1,3)
plot(ttm,sb3c,'b-',ttm,sb3,'r-')
datetick('x','mm/dd','keeplimits','keepticks');
ylabel('S (psu)','FontSize',14);
subplot(4,1,4)
plot(ttm,sb4c,'b-',ttm,sb4,'r-')
datetick('x','mm/dd','keeplimits','keepticks');
ylabel('S (psu)','FontSize',14);
xlabel('Date in 2009 (PST)','FontSize',14);

save saldata_bakerbay_comp2_32w.out saldata -ASCII -TABS

!rm hgrid2.gr3
!rm hgrid2.ll
!rm hgrid2c.gr3
!rm hgrid2c.ll
!rm *_????.6?

