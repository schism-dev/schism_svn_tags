function corie = corie_ywdn2corie(year, week, day, tsn, ntsDay)
% corie = corie_ywdn2corie(yyyy, ww, dd, tsn [, ntsDay])
% maps year, week number - ww, day of the week - dd, and timestep tsn
% to corie - corie date
% all in/out are simple numbers
% ntsDay number of timesteps a day (defualt is 96 (evry 15 mins))
%
% Sergey Frolov, November 2004

if nargin <5 
    ntsDay = 96;
end

%integer part of the day 
corie = datenum(['dec-31-',num2str(year-1)])-datenum('dec-31-1995')+((week-1)*7+(day));

%fractional part of the day
corie = corie + tsn/ntsDay;
