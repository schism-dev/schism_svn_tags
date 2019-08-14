function []=stVis(st,stateInfo,dbInfo,gr,flag,ln,maxmin) 
%[]=stVis(st,stateInfo,dbInfo,gr,flag,ln,[maxmin]) 
%visualise 2d slabs at level 'ln' from the state vector 'st'
%flag = ['e','s','t','u','v']
%ln = vertical layer number
%maxmin = [min max ] limits for color scale in units of st
%
%Sergey Frolov


if strcmp(flag,'e')
	dv=double(st(eval(stateInfo.idx.str.e),1))*stateInfo.var.e;
	sprintf('range is max: %g min: %g', max(dv),min(dv))
	if nargin == 7 
        gr_plot(gr.hgrid,dv,maxmin)
    else 
        gr_plot(gr.hgrid,dv)
    end
elseif strcmp(flag,'s')
	d3d.h   = dbInfo.fheaders.s;
	d3d.data= double(st(eval(stateInfo.idx.str.s),1));
	dv      = map_eb2hts(d3d)'*stateInfo.var.s;
	sprintf('range is max: %g min: %g', max(dv(:,ln)),min(dv(:,ln)))
	gr_plot(gr.hgrid,dv(:,ln),[min(dv(:)),max(dv(:))])
	if nargin == 7 
     	gr_plot(gr.hgrid,dv(:,ln),maxmin)
    else
        gr_plot(gr.hgrid,dv(:,ln))
    end
elseif strcmp(flag,'t')
	d3d.h   = dbInfo.fheaders.s;
	d3d.data= double(st(eval(stateInfo.idx.str.t),1));
	dv      = map_eb2hts(d3d)'*stateInfo.var.t;
	sprintf('range is max: %g min: %g', max(dv(:,ln)),min(dv(:,ln)))
	if nargin == 7 
     	gr_plot(gr.hgrid,dv(:,ln),maxmin)
    else
        gr_plot(gr.hgrid,dv(:,ln))
    end
elseif strcmp(flag,'u')
	d3d.h   = dbInfo.fheaders.s;
	d3d.data= double(st(eval(stateInfo.idx.str.u),1));
	dv      = map_eb2hts(d3d)'*stateInfo.var.u;
	sprintf('range is max: %g min: %g', max(dv(:,ln)),min(dv(:,ln)))
	if nargin == 7 
     	gr_plot(gr.hgrid,dv(:,ln),maxmin)
    else
        gr_plot(gr.hgrid,dv(:,ln))
    end
elseif strcmp(flag,'v')
	d3d.h   = dbInfo.fheaders.s;
	d3d.data= double(st(eval(stateInfo.idx.str.v),1));
	dv      = map_eb2hts(d3d)'*stateInfo.var.v;
	sprintf('range is max: %g min: %g', max(dv(:,ln)),min(dv(:,ln)))
	if nargin == 7 
     	gr_plot(gr.hgrid,dv(:,ln),maxmin)
    else
        gr_plot(gr.hgrid,dv(:,ln))
    end
else
    error('unknown flag')
end
    
