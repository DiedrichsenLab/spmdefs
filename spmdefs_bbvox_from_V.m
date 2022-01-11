function [bb,vx] = spmdefs_bbvox_from_V(V)
% Return the default bounding box for an image volume

vx = sqrt(sum(V.mat(1:3,1:3).^2));
o  = V.mat\[0 0 0 1]';
o  = o(1:3)';
bb = [-vx.*(o-1) ; vx.*(V.dim(1:3)-o)];
return;
