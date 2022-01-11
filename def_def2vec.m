function [VEC,V,D]=def_def2vec(name,varargin);
% function [VEC,V,D]=spmj_defs2vec(name,vecname,lengthname);
% Generates deformation vectors from deformation map
%   saves them as image files if names are given 
% vararginoptions:
%   - 'normalize',maskimg    : Make the mean over the displacement map within
%                         the mask zero, if mask empty over whole image
%  - 'affine', A : Subtract affine transformation from the vectors 
%       The affine matrix in the _sn.mat is from template to source voxel!
%       So, this need to be called with 
%       'affine',VF.mat*Affine*inv(VG.mat)
%  - 'vecname',name: Name of the vector output file 
%  - 'lengthname',name: Name of the length output file 
normalize=NaN;
affine=eye(4); 
vecname=[]; 
vararginoptions(varargin,{'normalize','affine','vecname'});

% Read in the deformation field 
for i=1:3
    V(i)=spm_vol(sprintf('%s,1,%d',name,i));
    D{i}=spm_read_vols(V(i));
end;

% read the normalisation volume 
if ~isnan(normalize)
    if ischar (normalize)
        VM=spm_vol(normalize);
        XM=spm_read_vols(VM);
    else 
        XM=ones(size(D{1}));
    end;
    XM=find(XM);
end;

% generate Vector map 
% By subtracting the voxel location 
% in the VG image from the location in the new image 
d=V(1).dim;
[I{1},I{2},I{3}]=meshgrid([1:d(1)],[1:d(2)],[1:d(3)]);
[I{1},I{2},I{3}]=spmj_affine_transform(I{1},I{2},I{3},affine*V(1).mat);
for i=1:3
    VEC{i}=D{i}-permute(I{i},[2 1 3]);
    if (~isnan(normalize)) 
        VEC{i}=VEC{i}-mean(VEC{i}(XM));
    end;
end;

if (~isempty(vecname))
    for i=1:3 
        Vvec(i)=V(i);
        Vvec(i).fname=vecname;
        Vvec(i)=spm_write_vol(Vvec(i),VEC{i});
    end;
end;