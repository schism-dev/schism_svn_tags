function [grid]=gr_compute_sc(grid,sccFname ,scFname )
%[grid]=gr_compute_sc(grid,sccFname [,scFname ])
% load indexes of side centers
% grid - existing grid structure as outputed by readGrid
% sccFname - side centers connectivity .bp file [ss_node nodeIdx1 nodeIdx2 junk]
% scFname - optional side centers .bp file [ss_node x y depth]
%
% Sergey Frolov March 08, 2004


if nargin == 3
	[scc]                   = gr_readGrid(sccFname);    
	[sc]                    = gr_readGrid(scFname);
	grid.sideCent.nodes     = sc.nodes;
	grid.sideCent.connect   = scc.nodes(:,1:3);
	grid.sideCent.np        = scc.np;
    grid.sideCent.sccFname  = scc.fname;
    grid.sideCent.scFname   = sc.fname;
    grid.sideCent.flag      = 1;
elseif nargin == 2
	[scc]                   = gr_readGrid(sccFname);    
	grid.sideCent.nodes     = [];
	grid.sideCent.connect   = scc.nodes(:,1:3);
	grid.sideCent.np        = scc.np;
    grid.sideCent.sccFname  = scc.fname;
    grid.sideCent.scFname   = [];
    grid.sideCent.flag      = 1;
end
    