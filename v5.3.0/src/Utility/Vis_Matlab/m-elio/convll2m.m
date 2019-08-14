function [x,y]=convll2m(lon,lat,RefLon,RefLat) 
%CONVLL2M convert long/lat (deg) to cartesian coordinates (meters)
% CONVLL2M converts from long/lat coordinates (in degrees) to
% meters.  
%
%    Input: lon  - longitudes in degrees
%           lat  - latitudes in degrees
%           RefLon - reference longitude (defaults to Boston, MA)
%           RefLat - reference latitude (defaults to Boston, MA)
%  Outputs: x   - east coordinate in meters
%           y   - west coordinate in meters
%
%  Call as:  [x,y]=convll2m(lon,lat,RefLon,RefLat);
%
% Written by : Brian O. Blanton
%              Summer 2000

x0 = lon;
y0 = lat;
lon = lon(:);
lat = lat(:);

if nargin == 3
    switch true
        case strcmp(RefLon,'Atlantic')
            method = 'blan';
            RefLon = 0;
            RefLat = 35.3256;
            R=6367500;  % Approx earth radius in meters
        case strcmp(RefLon,'Pacific')|nargin==2
            fid = fopen('temp.dat','w');
            if size(lon,1)==1
                data = [(1:length(lon));lon;lat;ones(1,length(lon))];
            else
                data = [(1:length(lon));lon';lat';ones(1,length(lon))];
            end 
            fid = fopen('temp.dat','w');
            fprintf(fid,'%s\n','***');
            fprintf(fid,'%i\n',size(lon,1));
            fprintf(fid,'%12i %12.3f %12.3f %i\n',data);
            fclose(fid);
            method = 'cmop';
    end
end
if nargin == 2
	if size(lon,1)==1
        data = [(1:length(lon));lon;lat;ones(1,length(lon))];
    else
        data = [(1:length(lon));lon';lat';ones(1,length(lon))];
    end 
	fid = fopen('temp.dat','w');
	fprintf(fid,'%s\n','***');
	fprintf(fid,'%i\n',size(lon,1));
	fprintf(fid,'%12i %12.4f %12.4f %i\n',data);
	fclose(fid);
	method = 'cmop';
end

if nargin < 2
   error('Too few arguments to CONVLL2M')
end
if nargout ~=2
   error('CONVLL2M requires two return arguments.')
end

% Check sizes of lon,lat
if ~all(size(lon)==size(lat))
   error('Sizes of lon and lat are NOT equal')
end

switch 1
    
    case strcmp(method,'blan')
        deg2rad=pi/180;
        fac=R*cos(RefLat*deg2rad);
        x=fac*lon*deg2rad;
        y=fac*log((1+sin(lat*deg2rad))./(cos(lat*deg2rad)));
        
    case strcmp(method,'cmop')
        !./spcs2ll_bp  -input   temp.dat    -output    temp2.dat    -ll2spcs
        fid = fopen('temp2.dat','r');
        out = fgets(fid);
        out = fgets(fid);
        out = fscanf(fid,'%i %f %f %f\n',[4 inf]);
        x = out(2,:);
        y = out(3,:);
        !rm -f temp.dat temp2.dat
        
end

x0(:) = x;
y0(:) = y;
x = x0;
y = y0;