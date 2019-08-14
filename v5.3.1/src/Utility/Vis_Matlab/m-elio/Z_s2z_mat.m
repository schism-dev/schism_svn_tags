function [z_rho,z_w] = Z_s2z_mat(h,zeta,S)
%
% Z_s2z_mat.m  6/12/2006  Parker MacCready
%
% This gives up to 3D arrays of z position (m, positive up)
% packed from the bottom to the top, for a given bottom depth (h)
% and surface height
% Note that the inputs "h" and "zeta" are matrices [LxM],
% and "S" is a structure
% created by "Z_get_basic_info.m" which contains all the s-coordinate info.
% z_rho is at box mid points, and z_w is at box top or bottom edges

% this code loops through vertical position

% get z positions from the S coordinates
% general version, loop through vertical levels
omat = ones(size(h));
for ii = 1:S.N
    zr0 = (S.s_rho(ii)-S.Cs_r(ii))*S.hc*omat + S.Cs_r(ii)*h;
    z_rho(ii,:,:) = zr0 + zeta.*(omat + zr0./h);
end
for ii = 1:S.N+1
    zw0 = (S.s_w(ii)-S.Cs_w(ii))*S.hc*omat + S.Cs_w(ii)*h;
    z_w(ii,:,:) = zw0 + zeta.*(omat + zw0./h);
end
