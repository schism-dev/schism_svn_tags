function [hts]=mapFile_eb2hts(gr, fname, itur, time, iths, ifile, eb_eta2, eb_snd, eb_ssd, eb_tnd, eb_tsd, eb_usg,tnd0,snd0)
%[hts]=mapFile_eb2hts(gr, fname, itur, time, iths, ifile, eb_eta2, eb_snd, eb_ssd, eb_tnd, eb_tsd, eb_usg,tnd0,snd0)
%map elcirc binary outout to hotstart files
%gr -grid as returned by readGrid
%data_np data at nodes data_ne data at elem
%
%Sergey Frolov march 2004

[hts]=hotstart_ini_da(gr,fname,itur);

hts.fid     = 0;
hts.time    = time;
hts.iths    = iths;

hts.eta1    = eb_eta2.data';
hts.eta2    = eb_eta2.data';
hts.we      = zeros(hts.nlev+1,hts.ne);

eb_usg_tmp  = eb_usg;
eb_usg_tmp.data     = eb_usg.data(:,1);
hts.ugsc            = map_eb2hts(eb_usg_tmp,0);
eb_usg_tmp.data     = eb_usg.data(:,2);
hts.vgsc            = map_eb2hts(eb_usg_tmp,0);

hts.tsd     = map_eb2hts(eb_tsd,1);
hts.ssd     = map_eb2hts(eb_ssd,1);
hts.tnd     = map_eb2hts(eb_tnd,1);
hts.snd     = map_eb2hts(eb_snd,1);

hts.tnd0    = tnd0;
hts.snd0    = snd0;

if itur
	hts.q2      = zeros(hts.nlev+1,hts.ns);
	hts.xl      = zeros(hts.nlev+1,hts.ns);
end
    
hts.ifile       = ifile;
hts.ifile_char  = sprintf('%12s',num2str(ifile));