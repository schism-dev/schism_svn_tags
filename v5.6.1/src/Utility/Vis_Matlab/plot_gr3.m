function []=plot_gr3(fname)
%plot_gr3('hgrid.gr3')
%Plot depths in .gr3 (pure tri) in matlab
%Adjuts color etc below

fid=fopen(fname,'r');
char=fgetl(fid);
tmp1=str2num(fgetl(fid));
ne=fix(tmp1(1));
np=fix(tmp1(2));
bathy(1:np,1:4)=nan;
nm(1:ne,1:5)=nan;
for i=1:np
  tmp=str2num(fgetl(fid));
  bathy(i,1:4)=tmp(:);
end %for i
for i=1:ne
  nm(i,:)=str2num(fgetl(fid));
end %for i
fclose(fid);

close all;
figure;
hold on;
%plot(bnd(:,1),bnd(:,2),'k.','MarkerSize',2);
patch('Faces',nm(:,3:5),'Vertices',bathy(:,2:3),'FaceVertexCData',bathy(:,4),'FaceColor','interp','EdgeColor','none');
colormap(jet(40));
%colormap(redblue);
caxis([-6 0]);
colorbar;
%axis([5e5 6e5 4.13e6 4.24e6]);
set(gcf,'Color',[1 1 1]);

