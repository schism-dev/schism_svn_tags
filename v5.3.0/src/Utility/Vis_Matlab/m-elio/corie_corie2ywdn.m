function [y, w, d, n] = corie_corie2ywdn(corie)
%[y, w, d, n] = corie_corie2ywdn(corie)
% converts corie date to 
% y-year, w-week, d-day of the week, n-timestep in day
%
% Sergey Frolov, June 2005


dn	= datenum('12/31/1995') + corie;		%number of days
dv	= datevec(dn);					%date vector
y	= dv(1);					%year
jd	= dn - datenum(['01/01/' num2str(y)]) + 1;	%julian day
pd	= jd - floor(jd);				%partial day
jd	= floor(jd);
w	= floor((jd-1)/7) + 1;
d	= jd - (w-1)*7;
%n	= floor(96*pd);
n      = round(96*pd);


%spetial case of singularity tsn =0 is actually stored in the previous day
if n ==0
  n = 96;
  d = d-1;
  if d == 0
    d=7;
    w=w-1;
    if w==0
      w=53;
      y=y-1;
      if mod(y,4)==0&(~mod(y,100)==0|mod(y,400)==0) %leap year
        d = 2;
      else
        d = 1;
      end
    end
  end
end
