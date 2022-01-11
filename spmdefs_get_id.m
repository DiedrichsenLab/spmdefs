function [Def,mat] = spmdefs_get_id(job)
% Get an identity transform based on an image volume.
N   = nifti(job.space{1});
d   = [size(N.dat),1];
d   = d(1:3);
mat = N.mat;
Def = cell(3,1);
[y1,y2,y3] = ndgrid(1:d(1),1:d(2),1:d(3));
Def{1} = single(y1*mat(1,1) + y2*mat(1,2) + y3*mat(1,3) + mat(1,4));
Def{2} = single(y1*mat(2,1) + y2*mat(2,2) + y3*mat(2,3) + mat(2,4));
Def{3} = single(y1*mat(3,1) + y2*mat(3,2) + y3*mat(3,3) + mat(3,4));
