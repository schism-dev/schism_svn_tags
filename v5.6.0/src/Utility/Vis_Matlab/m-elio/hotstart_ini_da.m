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

hts.length.head   = 8 + 4;
hts.length.elem   = ((2+hts.nlev+1)*hts.ne)*8;
hts.length.sides  = (hts.nlev*4+2)*hts.ns*8;
hts.length.points = (hts.nlev*(8+8+8+8))*hts.np;
if itur
	hts.length.itur  = 2*(hts.nlev+1)*hts.ns*8;
else
    hts.length.itur  = 0;
end

hts.startPos.head   = 0;
hts.startPos.elem   = hts.startPos.head + hts.length.head;
hts.startPos.sides  = hts.startPos.elem + hts.length.elem;
hts.startPos.points = hts.startPos.sides + hts.length.sides;
hts.startPos.itur   = hts.startPos.points + hts.length.points;
hts.startPos.tail   = hts.startPos.itur + hts.length.itur;

hts.fname   = fname;