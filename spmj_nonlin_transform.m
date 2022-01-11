function [X2,Y2,Z2] = spmj_nonlin_transform(prm,x,y,z)
% Transfers voxel location from atlas template to original space
% This saves the work of resampling the raw data into MNI space if you
% would like your data analyze in original space
Tr = prm.Tr;
[ax,ay,az]=spmj_affine_transform(x,y,z,inv(prm.VG.mat));

BX = spm_dctmtx(prm.VG(1).dim(1),size(Tr,1),ax-1);
BY = spm_dctmtx(prm.VG(1).dim(2),size(Tr,2),ay-1);
BZ = spm_dctmtx(prm.VG(1).dim(3),size(Tr,3),az-1);
% if flags.preserve,
% 	DX = spm_dctmtx(prm.VG(1).dim(1),size(Tr,1),x-1,'diff');
% 	DY = spm_dctmtx(prm.VG(1).dim(2),size(Tr,2),y-1,'diff');
% 	DZ = spm_dctmtx(prm.VG(1).dim(3),size(Tr,3),z-1,'diff');
% end;
% d  = [flags.interp*[1 1 1]' flags.wrap(:)];


% 	C = spm_bsplinc(V(i),d);

for j=1:length(z)                             % Cycle over all points 
    for k=1:size(BZ,2) 
        A(:,:,k)=BX(j,:)'*BY(j,:)*BZ(j,k);
    end;
    X1(j,1)=ax(j)+sum(sum(sum(Tr(:,:,:,1).*A)));
    Y1(j,1)=ay(j)+sum(sum(sum(Tr(:,:,:,2).*A)));
    Z1(j,1)=az(j)+sum(sum(sum(Tr(:,:,:,3).*A)));
    
end;
[X2,Y2,Z2]=spmj_affine_transform(X1,Y1,Z1,prm.VF.mat*prm.Affine);
% prm.VG.mat\prm.VF.mat*