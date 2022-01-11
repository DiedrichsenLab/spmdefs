function def_vtk2def(name,outname,target);
% function V=def_vtk2def(name,outname,target);
% Makes a vtk-file from the LDDMM algorithm into a y_* file for
% deformations
if (nargin<3)
    target=spm_get(1,'*.img','Choose Traget file');
end;

V_target=spm_vol(target); % Get information about target image
T=readVTK(name);
if (length(T)==0)
    error(sprintf('Could not find %s',name));
end;
V.fname=outname;
V.dim=size(T);
V.dim(4)=16;
V.pinfo=[1;0;0];
if (any(V.dim(1:3)~=V_target.dim(1:3)))
    error('target and vtk file do not have same dimensions');
end;
V.mat=eye(4); % V_target.mat;

% write our three volumes for the deformation map
for i=1:3
    Vo{i}=V;
    Vo{i}.n=i;
    Vo{i}=spm_create_vol(Vo{i});
    spm_write_vol(Vo{i},T(:,:,:,i));
    spm_close_vol(Vo{i});
end;
hdrname=[outname(1:end-3) 'hdr'];
[hdr,swap]=spm_read_hdr(hdrname);
hdr.hist.origin(1:3)=[1 1 1];
spmj_write_hdr(hdrname,hdr,swap);
