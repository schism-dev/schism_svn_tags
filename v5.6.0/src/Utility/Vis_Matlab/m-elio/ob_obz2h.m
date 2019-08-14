function [ob]=ob_obz2h(ob);
% [ob]=ob_obz2h(ob);
% computes a sparce observation matrix ob.z.H based on pidx and w in ob.z
% assumes that ob.z.flag == 2 (weights are already computed)
%
% Sergey Frolov, December 2004

if ob.z.flag <2 
    error('ob.z.flag<2, hence no pidx and w info computed yet')
end

i   = ob.z.pidx;
numObs  = size(i,1);
j   = repmat([1:numObs]',1,2);   %assumes that pidx [numObs, 2], where 2 are 2 brackiting levels
i   = i(:);
j   = j(:);
w   = ob.z.w(:);
%dimensions of matrix H=nxm  where
%m is dimension of a state to observe (num levels here)
%n is a number of measurements
m   = ob.gr.vgrid.nLevels;
n   = numObs;

H   = sparse(j,i,w,n,m);

ob.z.H = H;
ob.z.flag = 3;