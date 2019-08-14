function [ob]=ob_obxy2h(ob);
% [ob]=ob_obxy2h(ob);
% computes a sparce observation matrix ob.xy.H based on pidx and w in ob.xy
% assumes that ob.xy.flag == 2 (weights are already computed)
%
% Sergey Frolov, December 2004

if ob.xy.flag <2 
    error('ob.xy.flag<2, hence no pidx and w info computed yet')
end

i   = ob.xy.pidx;
numObs  = size(i,1);
j   = repmat([1:numObs]',1,3);   %assumes that pidx [numObs, 3], where 3 is 3 points of a triangle
i   = i(:);
j   = j(:);
w   = ob.xy.w(:);
%dimensions of matrix H=nxm  where
%m is dimension of a state to observe
%n is a number of measurements
m   = ob.gr.hgrid.np;
n   = numObs;

H   = sparse(j,i,w,n,m);

ob.xy.H = H;
ob.xy.flag = 3;