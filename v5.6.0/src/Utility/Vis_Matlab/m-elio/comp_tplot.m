clear;clear all;clc;close all;

base_out_dir ='/home/workspace/ccalmr45/choj/for_grant/barb_sal_maps/';
start_day = datenum(2009, 1, 1, 0, 0, 0);

curr_run_links_dir='/home/workspace/local0/forecasts';
modnames={'fca200nb','dev'};

yr=2009;imn=8;idy=19;
corie_day = datenum(1995, 12, 31, 0, 0, 0);
curr_day=datenum(yr, imn, idy, 0, 0, 0)+2;
doy=curr_day-start_day;
ne=720;ne2=710;
nnd=14;nst=109;

for imod=1:length(modnames)
   figure(10+imod)
    set(gcf,'Position',[1500 200 800 600],'Color','w')  
   for it=1:nnd
    model_day = sprintf('%d-%03d', yr, doy+it);
    eval(['!ln -sf ' curr_run_links_dir '/' modnames{imod} '/' model_day '/data/ogi02_CT1956.dat_CTD ogi02b.inp'])
    eval(['!ln -sf ' curr_run_links_dir '/' modnames{imod} '/' model_day '/data/ogi02_CT2946.dat_CTD ogi02s.inp'])
    eval(['!ln -sf ' curr_run_links_dir '/' modnames{imod} '/' model_day '/process/ogi02_elcirc_n001.dat ogi02m.inp'])
    load './ogi02s.inp';load './ogi02b.inp';load './ogi02m.inp';
    o2s=ogi02s;o2b=ogi02b;o2m=ogi02m;
    %-- Read Observed data S,T 
    tmp0=o2s(:,1);tmp1=o2s(:,3);tmp2=o2s(:,5);tmp0b=o2b(:,1);tmp3=o2b(:,3);tmp4=o2b(:,5);
    %-- Read Modeled data
    nlg=length(o2m);
    for kk=1:nlg/nst
        fnline=(kk-1)*nst;
        tmm(kk)=o2m(fnline+1,1);
        dp1s(kk)=o2m(fnline+53,2);dp6s(kk)=o2m(fnline+47,2);
        dp1t(kk)=o2m(fnline+53,3);dp6t(kk)=o2m(fnline+47,3);
    end
    
    subplot(211)
    title('Salinity Comparison at ogi02','FontSize',14)
    plot(tmp0+corie_day-1,tmp1,'b-',tmp0b+corie_day-1,tmp3,'r-');hold on;
    plot(tmm+corie_day-1,dp1s,'g-','LineWidth',2);hold on;
    plot(tmm+corie_day-1,dp6s,'k-','LineWidth',2);hold on;
    datetick('x','mm/dd');hold on;
    text(curr_day-4,32,modnames(imod),'FontSize',14)
%    set(gca,'XMinorTick','on')
%    xlim([curr_day curr_day+nnd+1])
%    set(gca,'XTick',[curr_day:2:curr_day+nnd+1])
%    set(gca,'XTickLabel',{'08/18','08/20','08/22','08/24','08/26','08/28','08/30','09/01','09/03'})
    grid on;
    legend('obs(1m)','obs(6m)','pre(1m)','pre(6m)',3)

    subplot(212)
    title('Temperature Comparison at ogi02','FontSize',14)
    plot(tmp0+corie_day-1,tmp2,'b-',tmp0b+corie_day-1,tmp4,'r-');hold on;
    plot(tmm+corie_day-1,dp1t,'g-','LineWidth',2);hold on;
    plot(tmm+corie_day-1,dp6t,'k-','LineWidth',2);hold on;
    datetick('x','mm/dd');hold on;
    text(curr_day-4,16,modnames(imod),'FontSize',14)
%    set(gca,'XMinorTick','on')
%    xlim([curr_day curr_day+nnd+1])
%    set(gca,'XTick',[curr_day:2:curr_day+nnd+1])
%    set(gca,'XTickLabel',{'08/18','08/20','08/22','08/24','08/26','08/28','08/30','09/1'})
    grid on;

%    !rm *.inp
   end
   if 0>1
    figure(imod)
    set(gcf,'Position',[1500 200 800 600],'Color','w')  
    subplot(211)
    plot(tmm+corie_day+2,sals,'b-',tmm+corie_day+2,salb,'r-')
    plot(tmm1+corie_day+2,sals1,'b-',tmm+corie_day+2,salb1,'r-')
    datetick('x','mm/dd')
    xlim([curr_day curr_day+16])
    set(gca,'XTick',[curr_day:2:curr_day+16])
    set(gca,'XTickLabel',{'08/16','08/18','08/20','08/22','08/24','08/26','08/28','08/30','09/1'})
    grid on;
    subplot(212)
    plot(tmm+corie_day+2,tems,'b-',tmm+corie_day+2,temb,'r-')
    datetick('x','mm/dd')
    xlim([curr_day curr_day+16])
    set(gca,'XTick',[curr_day:2:curr_day+16])
    set(gca,'XTickLabel',{'08/16','08/18','08/20','08/22','08/24','08/26','08/28','08/30','09/1'})
    grid on;
  end
end
