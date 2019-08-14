function [elAreas]=gr_compute_elAreas(gr)
%[elAreas]=gr_compute_elAreas(gr)

ec = gr.hgrid.elem(:,3:6);
idx3el = gr.hgrid.elem(:,2)==3;
idx4el = gr.hgrid.elem(:,2)==4;

elAreas=zeros(gr.hgrid.ne,1);
elAreas(idx3el,:)=polyarea(gr.hgrid.x(ec(idx3el,1:3))',gr.hgrid.y(ec(idx3el,1:3))');
elAreas(idx4el,:)=polyarea(gr.hgrid.x(ec(idx4el,:))',gr.hgrid.y(ec(idx4el,:))');

