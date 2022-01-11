function [x1,y1,z1] = spmdefs_transform(Def,mat,x,y,z,varargin)
% function [x1,y1,z1] = spmdefs_transform(Def,mat,x,y,z,varargin)
% Transforms coordinates from VG space into VF space 
% INPUT: 
%   Def: Deformation cell array {3} 
%   mat: affine transform for Def volume
%   x,y,z: source locations in world coordinates 
% VARAGRGIN
%   'interp',1: Default trialinear  

interp=1; 
vararginoptions(varargin,{'interp'}); 
[vx,vy,vz]=spmj_affine_transform(x,y,z,inv(mat)); 
x1=spm_sample_vol(Def{1},vx,vy,vz,interp); 
y1=spm_sample_vol(Def{2},vx,vy,vz,interp); 
z1=spm_sample_vol(Def{3},vx,vy,vz,interp); 
