function [hts]=hotstart_read_elcirc(fname,gr,itur)
%[hts]=hotstart_read_elcirc(fname,gr,itur)
% read elcirc hotstart file fname as defined in elcirc v5.02 
% gr - a computation grid as returned by readGrid
% itur[0|1] woth our without turbulence closure params
%
%Sergey Frolov March 2004

tic
%time,iths,(eta1(i),eta2(i), (we(i,j),j=0,nvrt),i=1,ne), ((vn2(i,j),vt2(i,j),j=0,nvrt),(tsd(i,j),ssd(i,j),j=1,nvrt),i=1,ns) ,(peta(i),ibad(i),(uu1(i,j),vv1(i,j),ww1(i,j),nosm(i,j),j=0,nvrt), 
%(tnd(i,j),snd(i,j),j=1,nvrt),i=1,np),((q2(i,j),xl(i,j),j=0,nvrt), i=1,ns),(irec(i),noutc(i),i=1,noutput),igmp,nscougm
% xl_side is optional (use if itur=3)

hts.ne      = gr.hgrid.ne;
hts.np      = gr.hgrid.np;
hts.nlev    = gr.vgrid.nLevels;
hts.ns      = gr.sideCent.np;
hts.itur    = itur;

hts.length.head   = 8 + 4;
hts.length.elem   = ((2+hts.nlev+1)*hts.ne)*8;
hts.length.sides  = (hts.nlev*4+2)*hts.ns*8;
hts.length.points = (8+4 + (hts.nlev+1)*(8+8+8+4) +hts.nlev*(8+8+8+8))*hts.np;
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
hts.fid     = fopen(fname,'r','l');

hts = readHead(hts);
hts = readElem(hts);
hts = readSides(hts);
hts = readPoints(hts);
hts = readItur(hts);
hts = readTail(hts);

hts.fid = fclose(hts.fid);
disp(['Read file ' fname ' at ' num2str(toc) 'secs'])

function [hts]=readHead(hts)
    fseek(hts.fid, hts.startPos.head, 'bof');
	hts.time    = fread(hts.fid,1,'double');
	hts.iths    = fread(hts.fid,1,'int32');

function [hts] = readElem(hts)
	fseek(hts.fid, hts.startPos.elem, 'bof');
	len         = (2+hts.nlev+1)*hts.ne; % number of values of eta1,eta2, we(i,j)
	block       = fread(hts.fid,len,'double'); % write array
	block       = reshape(block,3+hts.nlev,hts.ne);
	hts.eta1    = block(1,:);
	hts.eta2    = block(2,:);
	hts.we      = block(3:end,:);

function [hts] = readSides(hts)
    fseek(hts.fid, hts.startPos.sides, 'bof');
	len         = (hts.nlev*4+2)*hts.ns;
	block       = fread(hts.fid,len,'double'); % write values
	block       = reshape(block,2+hts.nlev*4,hts.ns);
	hts.vn2     = block(1:2:2*hts.nlev+2,:);
	hts.vt2     = block(2:2:2*hts.nlev+2,:);
	hts.tsd     = block((2*hts.nlev+3):2:end,:);
	hts.ssd     = block((2*hts.nlev+4):2:end,:);


function [hts] = readPoints(hts)
    np      = hts.np;
    nlev    = hts.nlev;
	npStart = hts.startPos.points;
    
	fseek(hts.fid,npStart,'bof');
	peta   = fread(hts.fid,np,'double',hts.length.points/hts.np-8)';
	
	fseek(hts.fid,npStart+8,'bof');
	ibad   = fread(hts.fid,np,'int32',hts.length.points/hts.np-4)';
	
	for i=1:nlev+1
		fseek(hts.fid,npStart+8+4+28*(i-1),'bof');    
		uu(i,:)   = fread(hts.fid,np,'double',hts.length.points/hts.np-8)';
	
		fseek(hts.fid,npStart+8+4+8+28*(i-1),'bof');    
		vv(i,:)   = fread(hts.fid,np,'double',hts.length.points/hts.np-8)';
        
		fseek(hts.fid,npStart+8+4+8+8+28*(i-1),'bof');    
		ww(i,:)   = fread(hts.fid,np,'double',hts.length.points/hts.np-8)';
        
		fseek(hts.fid,npStart+8+4+8+8+8+28*(i-1),'bof');    
		nosm(i,:)   = fread(hts.fid,np,'int32',hts.length.points/hts.np-8)';
	end
	
	fseek(hts.fid,npStart+8+4+(nlev+1)*(8+8+8+4),'bof');    
	tsnd    = fread(hts.fid,[4*nlev np],[num2str(4*nlev) '*double'],hts.length.points/hts.np-nlev*4*8);
	tnd     = tsnd(1:4:end,:);
	snd     = tsnd(2:4:end,:);
	tnd0    = tsnd(3:4:end,:);
	snd0    = tsnd(4:4:end,:);
	
	hts.peta    = peta;
	hts.ibad    = ibad;
	hts.uu      = uu;
	hts.vv      = vv;
	hts.ww      = ww;
	hts.nosm    = nosm;
	hts.tnd     = tnd;
	hts.snd     = snd;
	hts.tnd0     = tnd0;
	hts.snd0     = snd0;      

function hts=readItur(hts)
	fseek(hts.fid, hts.startPos.itur, 'bof');
    ns      = hts.ns;
    nlev    = hts.nlev;
	if hts.itur
		len         = 2*(nlev+1)*ns;
		xl_array    = fread(hts.fid,len,'double'); %
		xl_array    = reshape(xl_array,2*nlev+2,ns);
		q2          = xl_array(1:2:end,:);
		xl          = xl_array(2:2:end,:);
	else 
		q2=[];
		xl=[];
	end;
	hts.q2  = q2; 
	hts.xl  = xl;

function [hts]=readTail(hts)
    fseek(hts.fid, hts.startPos.tail, 'bof');
	hts.ifile       = fread(hts.fid,1,'int32');
	hts.ifile_char  = sprintf('%12s',fread(hts.fid,12,'uchar'));    
    
function [hts] = readPoints_slow(hts)
% version with explicit loops
    np      = hts.np;
    nlev    = hts.nlev;
	peta        = zeros(1,np);
	ibad        = peta;
	uu          = zeros(nlev+1,np);
	vv          = uu;
	ww          = uu;
	nosm        = uu;
	tsnd        = zeros(nlev*4,np);
	for i=1:np
       peta(i)  = fread(hts.fid,1,'double'); %peta
       ibad(i)  = fread(hts.fid,1,'integer*4'); %ibad
       for j=1:nlev+1
          uu(j,i)   = fread(hts.fid,1,'double'); %u
          vv(j,i)   = fread(hts.fid,1,'double'); %v
          ww(j,i)   = fread(hts.fid,1,'double'); %w
          nosm(j,i) = fread(hts.fid,1,'integer*4'); %nosm
      end;
       tsnd(:,i)    = fread(hts.fid,4*nlev,'double'); % ts at node
	end;
	tnd     = tsnd(1:4:end,:);
	snd     = tsnd(2:4:end,:);
	tnd0    = tsnd(3:4:end,:);
	snd0    = tsnd(4:4:end,:);
	
	hts.peta    = peta;
	hts.ibad    = ibad;
	hts.uu      = uu;
	hts.vv      = vv;
	hts.ww      = ww;
	hts.nosm    = nosm;
	hts.tnd     = tnd;
	hts.snd     = snd;
	hts.tnd0     = tnd0;
	hts.snd0     = snd0;    
    