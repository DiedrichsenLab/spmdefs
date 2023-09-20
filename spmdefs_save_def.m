function spmdefs_save_def(Def,mat,ofname)
% function spmdefs_save_def(Def,mat,ofname)
% Saves a deformation field as a nii-file 
% INPUT: 
%   Def: deformation field 
%   mat: transformation matrix 
%   ofname: outfilename 
% 

% taken from spm_defs 
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
% John Ashburner
% $Id: spm_defs.m 991 2007-11-02 14:01:28Z john $
if isempty(ofname), return; end;

dim   = [size(Def{1},1) size(Def{1},2) size(Def{1},3) 1 3]; % was 1 3 
dtype = 'FLOAT32-BE';
off   = 0;
scale = 1;
inter = 0;
dat   = file_array(ofname,dim,dtype,off,scale,inter);

N      = nifti;
N.dat  = dat;
N.mat  = mat;
N.mat0 = mat;
N.mat_intent  = 'Aligned';
N.mat0_intent = 'Aligned';
N.intent.code = 'VECTOR';
N.intent.name = 'Mapping';
N.descrip = 'Deformation field';
create(N);
N.dat(:,:,:,1,1) = Def{1}; % was 1,1 
N.dat(:,:,:,1,2) = Def{2}; % was 1,2 
N.dat(:,:,:,1,3) = Def{3}; % was 1,3 
return;
