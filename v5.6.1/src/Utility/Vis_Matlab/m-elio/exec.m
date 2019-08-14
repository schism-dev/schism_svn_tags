function []=exec(cmd,echoflg,errorFlg)
% []=exec(cmd,[echoflg=0],[errorFlg=1])
% executes any shell cmd
% if errorFlg =1 fatal errors in cmd will raise error in matlab
% if errorFlg =0 and cmd fails cmd output is printed on the screen no error is raised
%    in any case output output of failed cmd is printed on the screen
% if echoflg =1 output of the program is prented on the screen
%
% sergey frolov March 11, 2005

if nargin <3
  errorFlg=1;
end
if nargin <2
  echoflg=0;
end


if echoflg==0
  [s,w]=system(cmd);
elseif echoflg==1
  [s,w]=system(cmd,'-echo');
end

if s~=0&errorFlg==1
   disp(w)
   error(['execution of ' cmd ' failed'])
elseif s~=0&errorFlg==0
   disp(w)
   warning(['execution of ' cmd ' failed'])
end

