function []=SELFE_TRANSECT(base,hgrid,transect_bp,binary_name,stacks,nspool,test)
% Authors: Sergey Frolov, Paul Turner, Joseph Zhang
% Date: April 2011
% Matlab function to visualize vertical transects
% If there is only 1 point in the transect input, vertical 'sample' is visualized
% If there are more than 1 point in the transect input, transect animaton is produced,
% and in this case for 3D vectors, currently only the x-component is plotted
% 
% Before calling this function, add path of m-elio library; e.g.
% path(path,'/usr/local/cmop/matlab/cmop/m-elio');
%
% SELFE_TRANSECT(base,hgrid,transect_bp,binary_name,stacks,nspool,test)
% Inputs: 
%         base: base directory where 'hgrid' and 'transect_bp' (below) reside (binary files are in base/outputs/)
%         hgrid: 'hgrid.gr3' or 'hgrid.ll';
%         transect_bp: build point file name defining the transect (see 'transect.bp.sample');
%         binary_name = 'salt.63' etc (no stack #)
%         stacks: array of stack numbers (e.g. [2 4 5]) in the output 
%                 file names (related to time)
%         nspool: sub-sampling frequency within each stack (e.g. 1 - include all)
%         test: 'y' (only plot out 1st frame for test); 'n' (plot all frames)
%         In addition, *elev.61 must also be inside base/outputs/
%         May need to adjust some parameters inside (e.g. caxis) to get right appearance of images
% Outputs: images and transect.avi (for transect only)

% Add m-elio to the matlab path (do this before calling any functions)
%path(path,'C:\yinglong\m-elio\m-elio');
%path(path,'~/bin/m-elio/m-elio/');
addpath 'C:\yinglong\m-elio' 'C:\yinglong\m-elio\m-elio';
addpath '~/bin/m-elio/' '~/bin/m-elio/m-elio/';

close all; 
scrsz = get(0,'ScreenSize'); %screen size
%figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]); 
figure('Position',[1 scrsz(4)*0.2 scrsz(3)/2 scrsz(4)*0.7]);

% Read a lat/lon .ll grid for other coordinates
tr.hgrid=gr_readHGrid(strcat(base,'/',transect_bp));
gr.hgrid=gr_readHGrid(strcat(base,'/',hgrid));
% Read the header for variable and vertical grid
header1=sz_readHeader(strcat(base,'/outputs/',num2str(stacks(1)),'_',binary_name));
gr.vgrid=header1.vgrid;
headerz=sz_readHeader(strcat(base,'/outputs/',num2str(stacks(1)),'_elev.61'));

dtout=header1.dt; %time step
nrec=header1.nSteps; % # of records in each stack
ivs=header1.flagSv; %1; scalar; 2: vector
i23D=header1.flagDm; %2 or 3D
nvrt=header1.vgrid.nLevels; %surface level index

% compute transect
[ob]= ob_ini_fromTrasect(gr,strcat(base,'/',transect_bp));
trLen = cumsum(sqrt((ob.xy.x(2:end)-ob.xy.x(1:end-1)).^2+(ob.xy.y(2:end)-ob.xy.y(1:end-1)).^2));
trLen = [0; trLen]; %trLen(1:nbp) is x-coord. of the transect
nbp=length(trLen);
% Transform depth vector to along-transect depths
dp=ob.xy.H*gr.hgrid.depth; 

% Sanity checks
if(i23D==2)
  error('2D var has no transect');
end
if(ivs==1 && nbp==1 && strcmp(test,'y'))
  error('Please set test to n and retry');
end

% Axis limits
xmax=max(trLen);
xmin=min(trLen);
ymax=0;
ymin=-30;

% Read the variable and the vertical grid
delete('transect.avi');
%avi_out = avifile('transect.avi','FPS',5); %output movie 
vidObj = VideoWriter('transect.avi');
open(vidObj);

if(strcmp(test,'y'))
  stacks2=stacks(1); it2=1;
else
  stacks2=stacks; it2=nrec;
end

time_min=(stacks2(1)-1)*nrec*dtout/3600; %for axis
time_max=stacks2(end)*nrec*dtout/3600;
count=0;
for day=stacks2
  timeout0=(day-1)*nrec*dtout;

  %Read binary files
  header1=sz_readHeader(strcat(base,'/outputs/',num2str(day),'_',binary_name));
  headerz=sz_readHeader(strcat(base,'/outputs/',num2str(day),'_elev.61'));

  %Outputs from the function sz_readTimeStep:
  %ts: info for each record (time stamp, iteration #, and elevation);
  %d(dataVector,vectorComponent_k,time_iteration)
  % where: dataVector is of size 'gridSize' (sum(np*(local levels)) for 3D);
  %        vectorComponent_k is either 1 or 2 (scalar or vector);
  %        time_iteration is record number in the binary file.
  [d ts]=sz_readTimeStep(header1,1:nrec); 
  [dz tsz]=sz_readTimeStep(headerz,1:nrec); 

  for it=1:nspool:it2;
    count=count+1;
    % Time info
    timeout=timeout0+it*dtout;
    time_d=fix(timeout/86400);
    time_h=fix((timeout-time_d*86400)/3600);
    time_m=fix((timeout-time_d*86400-time_h*3600)/60);
    time_s=timeout-time_d*86400-time_h*3600-time_m*60;

    % For 3D variables, map outputs to a regular matlab structure
    % Last argument "1" is a flag for vertical levels 1:nvrt (0 for ELCIRC).
    % Output u is a simple 2d aray (1:np,1:nvrt) for 3D variables
    [u] = map_sz2hts(header1,d(:,1,it),1);
    % Transform u onto along-transect array: uout(1:nbp,1:nvrt)
    uout=ob.xy.H*double(u); 

    if(ivs==2)
      [v] = map_sz2hts(header1,d(:,2,it),1);
      vout=ob.xy.H*double(v);
    end

    % Construct vertical grid
    % Read timestep of elevations
    eta=ob.xy.H*double(dz(:,1,it));
    sz=sz_computeZlevels(dp,eta,gr.vgrid); %sz(1:nbp,1:nvrt) is the y-coord. for the transect

    % plot
    if(nbp==1) %sample
      %Contruct matrix for scalars (plot later; inside it loop to conserve memory)
      sout(:,count)=uout(1,:);
      if(ivs==2); tout(:,count)=vout(1,:); end;
      sxout(:,count)=timeout/3600*ones(nvrt,1); %x coord.
      syout(:,count)=sz(1,:); %y coord.
      continue;
    else %transect
      hold on;
      axis([xmin xmax ymin ymax]);
      v2=axis;
      % Write time stamp info
      loc_info_x=(v2(2)+v2(1))/2;
      loc_info_y=v2(4)*1.05-v2(3)*0.05;
      text(loc_info_x,loc_info_y,{'Time (DD:HH:MM:SS)'; num2str([time_d time_h time_m time_s])});

      h=pcolor(trLen, sz', uout');
      set(h,'EdgeColor','none','FaceColor','interp');
      colormap(jet(15));
      caxis([0 30]); colorbar;
      xlabel('Along transect (m)');
      ylabel('z (m)');
      % Add image to avi file
%      frame = getframe(gcf);
%      avi_out=addframe(avi_out,frame);

      set(gca,'nextplot','replacechildren');
      currFrame = getframe(gcf);
      writeVideo(vidObj,currFrame);
      if(day ~= stacks2(end) || it ~=it2)
        clf; %clear figure to avoid overlay
      end
    end %nbp
  end %it
end %for day=stacks2

%Finish plotting for sample
%sout(nvrt,count) (count is # of time steps); tout (for vectors)
if(nbp==1)
  axis([time_min*0.99 time_max*1.01 ymin ymax]);
  hold on;

  if(ivs==1) 
    h=pcolor(sxout,syout,sout);
    set(h,'EdgeColor','none','FaceColor','interp');
    colormap(colorcube(40));
    caxis([0 20]); colorbar;
  else %vector
    quiver([sxout; sxout(end,:)],[syout; (ymin*0.98+ymax*0.02)*ones(1,count)], ...
    [sout; 1*ones(1,count)],[tout; 0*ones(1,count)],0,'k');
    text(time_min*0.7+time_max*0.3,ymin*0.98+ymax*0.02,'1 m/s');
  end %ivs
  xlabel('Time (hours)');
  ylabel('z (m)');
end %sample plot

%avi_out=close(avi_out);
close(vidObj);

