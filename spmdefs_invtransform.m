function [x1,y1,z1] = spmdefs_invtransform(coord,deformation,varargin)
% transforms subjects xyz locations into atlas space  
% INPUT: 
%   coord: Px3 matrix for x,y,z coordinates in individual subject space 
%   deformation: filename of the subj->atlas deformation 
% VARARGIN: 
%   'interp',0: interpolation (default 0: nearest neighbor) 
% OUTPUT: 
%   atlasCoord: Px3 matrix for x,y,z coordinates in atlas space 
intrp = 0;
vararginoptions(varargin,{'intrp'});
%_________________________________________________________________________
% get the deformation map in x y z. tells you the position in the original
% image as a function of position in suit space
[def,def_mat] = spmdefs_get_sn2def(deformation);

T=load(deformation);
VSource=T.VF;
[i_def,i_defMat] = spmdefs_get_inv(def,def_mat,VSource);
%voxelspace of the anatomical!!!!
[i j k]= spmj_affine_transform(coord(:,1),coord(:,2),coord(:,3),inv(VSource.mat));
for n=1:3 
    atlasCoord(:,n)=spm_sample_vol(i_def{n},i,j,k,interp);  
end; 
