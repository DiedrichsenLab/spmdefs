function [Def,mat] = spmdefs_get_inv(Def0,mat0,VT)
% function [Def,mat] = spmdefs_get_inv(Def0,mat0,VT)
% Invert a deformation field (derived from a composition of deformations)
M0      = mat0;
M1      = inv(VT.mat);
M0(4,:) = [0 0 0 1];
M1(4,:) = [0 0 0 1];
[Def{1},Def{2},Def{3}]    = spm_invdef(Def0{:},VT.dim(1:3),M1,M0);
mat         = VT.mat;
