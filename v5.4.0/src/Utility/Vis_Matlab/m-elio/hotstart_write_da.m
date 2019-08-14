function []=hotstart_write_da(fname,hts)
%[]=hotstart_write(fname,hts)
% write elcirc hotstart structure hts into a file fname
%
%Sergey Frolov March 2004

tic

hts.fid     = fopen(fname,'w','l');

writeHead(hts);
writeElem(hts);
writeSides(hts);
writePoints(hts);
if hts.itur
	writeItur(hts);
end
writeTail(hts);

hts.fid = fclose(hts.fid);
%disp(['Wrote file ' fname ' in ' num2str(toc) 'secs'])

function writeHead(hts)
    fseek(hts.fid, hts.startPos.head, 'bof');
	c       = fwrite(hts.fid,hts.time,'double');
	c       = fwrite(hts.fid,hts.iths,'int32');

function writeElem(hts)
	fseek(hts.fid, hts.startPos.elem, 'bof');
    block   = [hts.eta1; hts.eta2; hts.we];
	c       = fwrite(hts.fid,block,'double');

function writeSides(hts)
    nlev                        = hts.nlev;
    block1                      = zeros([2*nlev+2,hts.ns]);
    block2                      = zeros([2*nlev,hts.ns]);    
    block1(1:2:2*(nlev+1),:)    = hts.vn2;
    block1(2:2:(2*(nlev+1)+1),:)= hts.vt2;    
    block2(1:2:2*(nlev),:)      = hts.tsd;
    block2(2:2:(2*nlev+1),:)    = hts.ssd;    
    block                       = [block1; block2];    
    fseek(hts.fid, hts.startPos.sides, 'bof');
	c           = fwrite(hts.fid,block,'double'); 

function writePoints(hts)
    np      = hts.np;
    nlev    = hts.nlev;
	npStart = hts.startPos.points;
    
    block                   = zeros([4*nlev,hts.np]);
    block(1:4:4*nlev,:)     = hts.tnd;
    block(2:4:(4*nlev+1),:) = hts.snd;    
    block(3:4:(4*nlev+2),:) = hts.tnd0;
    block(4:4:(4*nlev+3),:) = hts.snd0;
    
	fseek(hts.fid, hts.startPos.points, 'bof');    
	c       = fwrite(hts.fid,block,'double');

	
function writeItur(hts)
    block   = zeros([2*(hts.nlev+1),hts.ns]);
    block(1:2:2*(hts.nlev+1),:)     = hts.q2;
    block(2:2:(2*(hts.nlev+1)+1),:) = hts.xl;    
    fseek(hts.fid, hts.startPos.itur, 'bof');
    c       = fwrite(hts.fid,block,'double'); 

function writeTail(hts)
    fseek(hts.fid, hts.startPos.tail, 'bof');
	c   = fwrite(hts.fid,hts.ifile,'int32');
	c   = fwrite(hts.fid,hts.ifile_char,'char');    
    
    
%time,iths,(eta1(i),eta2(i), (we(i,j),j=0,nvrt),i=1,ne), ((vn2(i,j),vt2(i,j),j=0,nvrt),(tsd(i,j),ssd(i,j),j=1,nvrt),i=1,ns) ,(peta(i),ibad(i),(uu1(i,j),vv1(i,j),ww1(i,j),nosm(i,j),j=0,nvrt), 
%(tnd(i,j),snd(i,j),j=1,nvrt),i=1,np),((q2(i,j),xl(i,j),j=0,nvrt), i=1,ns),(irec(i),noutc(i),i=1,noutput),igmp,nscougm
% xl_side is optional (use if itur=3)
    
