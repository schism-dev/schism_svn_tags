function [sideCent] = gr_readSCgrid(sccFname,scFname,snxyFname)
%[grid] = gr_readSCgrid(sccFname,scFname,[snxyFname])
% read connectivity table for side centers
% grid - existing grid structure as outputed by readGrid
% sccFname - side centers connectivity .bp file [ss_node nodeIdx1 nodeIdx2 junk]
% scFname - side centers .bp file [ss_node x y dps]
% snxyFname - file with the dump of snx sny varaibles from ELCIRC, file is in .bp format
%
% Sergey Frolov, March 08, 2004
% Sergey Frolov, Feb 22, 2005, added snxy


[scc]               = gr_readHGrid(sccFname);    
[sc]                = gr_readHGrid(scFname);    
sideCent.nodes      = sc.nodes(:,1:4);
sideCent.connect    = scc.nodes(:,1:3);
sideCent.np         = scc.np;
sideCent.sccFname   = scc.fname;
sideCent.scFname    = sc.fname;
sideCent.ne         = 0;
sideCent.elem       = [];
sideCent.eofLines   = {};
sideCent.nn         = sc.nodes(:,1);
sideCent.x          = sc.nodes(:,2);
sideCent.y          = sc.nodes(:,3);
sideCent.depth      = sc.nodes(:,4);

if nargin == 3
    [snxy]          = gr_readHGrid(snxyFname);
    sideCent.snx    = snxy.nodes(:,2);
    sideCent.sny    = snxy.nodes(:,3);
end

sideCent.flag      = 1;

    
