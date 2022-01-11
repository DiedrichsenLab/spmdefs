function [dJ,mat]=spmdefs_jacobiandet(filename,varargin)
% function [dJ,mat]=spmdefs_jacobiandet(filename,varargin)
% INPUT: 
%   filename: deformation field / filename 
% OUTPUT: 
%   
mask=''; 
outfilename=''; 
vararginoptions(varargin,{'mask','outfilename'}); 

[dir,name,ext]=spm_fileparts(filename); 
if (ext=='.mat')
    [Def,mat]=spmdefs_get_sn2def(filename); 
    spmdefs_save_def(Def,mat,'temp');
    filename='y_temp.nii'; 
    ext='.nii';
end; 
STEPSIZE=5000; 

V=spm_vol(filename); 

[X,Y,Z]=ndgrid([1:V(1).dim(1)],[1:V(1).dim(2)],[1:V(1).dim(3)]); 
if (~isempty(mask)) 
    M=spm_vol(mask); 
    Mx=spm_read_vols(M); 
    I=find(Mx>0)'; 
    N=length(I);
else 
    N=size(X(:),1); 
    I=[1:N]; 
end; 

v=zeros(STEPSIZE,1);
d=zeros(3,3,STEPSIZE);
D=zeros(size(X(:),1),1)*NaN;
for i=1:3 
    V(i)=spm_vol(sprintf('%s,1,%d',filename,i)); 
end; 

K=ceil(N/STEPSIZE); 
for n=1:K; 
    indx=[(n-1)*STEPSIZE+1:min(N,n*STEPSIZE)]; 
    if (n==K)
        clear v d;
    end; 
    for i=1:3 
        [v,d(i,1,:),d(i,2,:),d(i,3,:)]=spm_sample_vol(V(i),X(I(indx)),Y(I(indx)),Z(I(indx)),1);
    end; 
    for j=1:length(indx); 
    D(I(indx(j)))=det(d(:,:,j)); 
    end; 
    fprintf('.'); 
end;
D=reshape(D,V(1).dim);
fprintf('\n'); 

% Write the image to disk 
if (isempty(outfilename))
    outfilename=fullfile(dir,[name '_detJ' ext]); 
end; 
VO=V(1); 
VO.fname=outfilename; 
VO.dt=[16 1]; 
spm_write_vol(VO,D); 