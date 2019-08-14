function [hts]=hotstart_read_da(fname,gr,itur)
%[hts]=hotstart_read_da(fname,gr,itur)
% read elcirc hotstart file fname !!!modified!!! format 
% gr - a computation grid as returned by readGrid
% itur[0|1] woth our without turbulence closure params
%
%Sergey Frolov March 2004

tic
%time,iths,(eta1(i),eta2(i), (we(i,j),j=0,nvrt),i=1,ne), ((da_ugsc(i,j),da_ugsc(i,j),j=0,nvrt),(tsd(i,j),ssd(i,j),j=1,nvrt),i=1,ns) ,
%((tnd(i,j),snd(i,j),j=1,nvrt),i=1,np),((q2(i,j),xl(i,j),j=0,nvrt), i=1,ns),(irec(i),noutc(i),i=1,noutput),ifile,ifile_char


hts     = hotstart_ini_local(fname,gr,itur);
hts.fid = fopen(hts.fname,'r','l');
	hts = readHead(hts);
	hts = readElem(hts);
	hts = readSides(hts);
	hts = readPoints(hts);
	hts = readItur(hts);
	hts = readTail(hts);
hts.fid = fclose(hts.fid);
disp(['Read file ' fname ' at ' num2str(toc) 'secs'])

function [hts]=hotstart_ini_local(fname,gr,itur)
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
	hts.ugsc    = block(1:2:2*hts.nlev+2,:);
	hts.vgsc    = block(2:2:2*hts.nlev+2,:);
	hts.tsd     = block((2*hts.nlev+3):2:end,:);
	hts.ssd     = block((2*hts.nlev+4):2:end,:);


function [hts] = readPoints(hts)
    np          = hts.np;
    nlev        = hts.nlev;
	npStart     = hts.startPos.points;
    	
	fseek(hts.fid,hts.startPos.points,'bof');    
	tsnd        = fread(hts.fid,[4*nlev np],[num2str(4*nlev) '*double'],hts.length.points/hts.np-nlev*4*8);
	hts.tnd     = tsnd(1:4:end,:);
	hts.snd     = tsnd(2:4:end,:);
	hts.tnd0    = tsnd(3:4:end,:);
	hts.snd0    = tsnd(4:4:end,:);

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
    
