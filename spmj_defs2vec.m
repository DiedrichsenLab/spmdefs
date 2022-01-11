function [VEC,LEN]=spmj_defs2vec(name,vecname,lengthname,varargin);
% function [VEC,LEN]=spmj_defs2vec(name,vecname,lengthname);
% Generates deformation vectors from deformation map
%   saves them as image files if names are given 
% vararginoptions:
%   - 'normalize',maskimg    : Make the mean over the displacement map within
%                         the mask zero, if mask empty over whole image
normalize=NaN;
vararginoptions(varargin,'normalize');

% Read in the deformation field 
for i=1:3
    V(i)=spm_vol(sprintf('%s,%d',name,i));
    D{i}=spm_read_vols(V(i));
end;

if ~isnan(normalize)
    if ischar (normalize)
        VM=spm_vol(normalize);
        XM=spm_read_vols(VM);
    else 
        XM=ones(size(D{1}));
    end;
    XM=find(XM);
end;

d=V(1).dim;
[I{1},I{2},I{3}]=meshgrid([1:d(1)],[1:d(2)],[1:d(3)]);
[I{1},I{2},I{3}]=spmj_affine_transform(I{1},I{2},I{3},V(1).mat);
for i=1:3
    VEC{i}=D{i}-permute(I{i},[2 1 3]);
    if (~isnan(normalize)) 
        VEC{i}=VEC{i}-mean(VEC{i}(XM));
    end;
    if (~isempty(vecname))
        Vvec(i)=V(i);
        Vvec(i).fname=vecname;
        spm_write_vol(Vvec(i),VEC{i});
    end;
end;
LEN=sqrt(VEC{1}.^2+VEC{2}.^2+VEC{3}.^2);
if (~isempty(lengthname))
    Vlen=V(1);
    Vlen.fname=lengthname;
    spm_write_vol(Vlen,LEN);
end;