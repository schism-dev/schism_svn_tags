function [hts]=hotstart_ini_da(gr,fname,itur)
%[hts]=hotstart_ini(gr,fname,otur)
%initialise a hotstart file
%itur == 1 with turbulance closure par
%
%sergey Frolov  March, 2004


hts.ne      = gr.hgrid.ne;
hts.np      = gr.hgrid.np;
hts.nlev    = gr.vgrid.nLevels;
hts.ns      = gr.sideCent.np;
hts.itur    = itur;

hts.startPos.head   = 0;
hts.startPos.elem   = hts.startPos.head + 8 + 4;
hts.startPos.sides  = hts.startPos.elem + ((2+hts.nlev+1)*hts.ne)*8;
hts.startPos.points = hts.startPos.sides + (hts.nlev*4+2)*hts.ns*8;
hts.startPos.itur   = hts.startPos.points + (8+4 + (hts.nlev+1)*(8+8+8+4) +hts.nlev*(8+8+8+8))*hts.np;
if itur
	hts.startPos.tail   = hts.startPos.itur + 2*(hts.nlev+1)*hts.ns*8;
else
    hts.startPos.tail   = hts.startPos.itur;
end

hts.fname   = fname;