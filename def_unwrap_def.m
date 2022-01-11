function def_unwrap_def(name);
% Unwraps a deformation file from LDDMM
if (nargin<1 | isempty(name))
    name=spm_get(3,'*.img','Get deformation map');
    V=spm_vol(name)
else
    fname=name;
    for i=1:3
        n{i}=sprintf('%s,%d',fname,i);
        V(i)=spm_vol(n{i});
    end;
    name=char(n);
end;
for i=1:3
    D{i}=spm_read_vols(V(i));
end;
d=V(1).dim;
[I{1},I{2},I{3}]=meshgrid([1:d(1)],[1:d(2)],[1:d(3)]);
[I{1},I{2},I{3}]=spmj_affine_transform(I{1},I{2},I{3},V(1).mat);

for i=1:3
    I{i}=permute(I{i},[2 1 3]);
    flipped_pos=find(D{i}-I{i}>d(i)/2);
    flipped_neg=find(D{i}-I{i}<-d(i)/2);
    D{i}(flipped_pos)=D{i}(flipped_pos)-d(i);
    D{i}(flipped_neg)=D{i}(flipped_neg)+d(i);
end;
for i=1:3
    Vo(i)=V(i);
    [a,b]=fileparts(name(i,:));
    Vo(i).fname=sprintf('%s/unw_%s.img',a,b);
    Vo(i)=spm_create_vol(Vo(i));
    spm_write_vol(Vo(i),D{i});
    Vo(i)=spm_close_vol(Vo(i));
end;
