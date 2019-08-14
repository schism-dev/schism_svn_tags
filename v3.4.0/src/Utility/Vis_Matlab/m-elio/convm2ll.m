function [lon,lat]=convm2ll(x,y,RefLon,RefLat);
%CONVM2LL convert cartesian coordinates (meters) to long/lat (deg)
% CONVM2LL converts from cartesian coordinates in meters to
% long/lat coordinates by inverting the conversion formulas 
% used in the routine CONVLL2M.  This is the inverse of
% converting long/lat coordinates to cartesian coordinates in
% meters. The x,y  coordinates MUST have been computed using
% CONVLL2M, and the  reference long/lat pair MUST be the same
% as those used in the  "forward" conversion (defaults to Boston, MA).
% Any translation in the x,y coordinates that might have been
% done MUST be "undone" prior to conversion back to long/lat
% coordinates.  The inversion for latitude is done with 5 
% Newton-Raphson iterations with an initial guess of the
% reference latitude.  The longitude is trivial.
%
%  Inputs: x,y - cartesian coordinates in meters
%          RefLon,RefLat - reference long/lat values
%
% Outputs: lon,lat - long/lat coordinates.
%
% Call as:  [lon,lat]=convm2ll(x,y,RefLon,RefLat);

% Implement default lat & lon

lon0 = x;
lat0 = y;
x = x(:);
y = y(:);

if nargin == 3
    switch 1
        case strcmp(RefLon,'Atlantic')
            method = 'blan';
            RefLon = 0;
            RefLat = 35.3256;
            R=6367500;  % Approx earth radius in meters
        case strcmp(RefLon,'Pacific')
            fid = fopen('temp.dat','w');
            if size(x,1)==1
                data = [(1:length(x));x;y;ones(size(y))];
            else
                data = [(1:length(x))';x';y';ones(size(x'))];
            end 
            fid = fopen('temp.dat','w');
            fprintf(fid,'%s\n','***');
            fprintf(fid,'%i\n',length(x));
            fprintf(fid,'%12i %12.3f %12.3f %i\n',data);
            fclose(fid);
            method = 'cmop';
    end
end
if nargin == 2
	if size(x,1)==1
        data = [(1:length(x));x;y;ones(size(y))];
    else
        data = [(1:length(x));x';y';ones(size(x'))];
    end 
	fid = fopen('temp.dat','w');
	fprintf(fid,'%s\n','***');
	fprintf(fid,'%i\n',length(x));
	fprintf(fid,'%12i %12.3f %12.3f %i\n',data);
	fclose(fid);
	method = 'cmop';
end

if nargin < 2
   error('Too few arguments to CONVM2LL')
end

if nargout ~=2
   error('CONVM2LL requires two return arguments.')
end

% Check sizes of lon,lat
if ~all(size(x)==size(y))
   error('Sizes of lon and lat are NOT equal')
end

switch 1
    
    case strcmp(method,'blan')
        deg2rad=pi/180;
        fac=R*cos(RefLat*deg2rad);

        % Invert for lon (easy)
        lon=x/fac/deg2rad;
        phin=RefLat.*ones(size(x)).*deg2rad;
        CC=exp(y/fac);
        
        % 5 Newton-Raphson iterations, assume convergence!!
        for i=1:5
            phinm1=phin;
            numer=CC - (1+sin(phinm1))./(cos(phinm1));
            denom= -1- (1+sin(phinm1)).*tan(phinm1)./cos(phinm1);
            phin=phinm1-numer./denom;
        end

        lat=phin/deg2rad;
        
    case strcmp(method,'cmop')
        !./spcs2ll_bp  -input   temp.dat    -output    temp2.dat    -spcs2ll
        fid = fopen('temp2.dat','r');
        out = fgets(fid);
        out = fgets(fid);
        data = fscanf(fid,'%i %f %f %f',[4 inf]);
        lon = data(2,:);
        lat = data(3,:);
        !rm -f temp.dat temp2.dat
   
end

lon0(:) = lon;
lat0(:) = lat;
lon = lon0;
lat = lat0;