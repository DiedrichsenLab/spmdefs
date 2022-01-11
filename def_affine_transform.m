function def_affine_transform(P,A,B,outfile);
% Does a affine transform on a deformation field
% deformation field z=f(x)
% new deformation field is z = A f(Bx)
% 
%______________________________________________________________________
% @(#)def_affine_transform.m	1.0 Joern Diedrichsen 24/7/2006

V   = spm_vol([repmat([P ','],3,1) num2str([1 2 3]')]);
hld = 1;

VO = V;
for i=1:length(VO),
	VO(i).fname = [outfile ',' num2str(i)];
	VO(i).desc  = 'Affine trans deformation field';
    VO(i).mat=VO(i).mat*inv(B);
end;

VO = spm_create_vol(VO);

spm_progress_bar('Init',VO(1).dim(3),'Transforming deformation','planes completed');

for p=1:VO(1).dim(3),
    % Slice first image
	M  = spm_matrix([0 0 p]);
	y1 = spm_slice_vol(V(1), M, V(1).dim(1:2),1);
	y2 = spm_slice_vol(V(2), M, V(1).dim(1:2),1);
	y3 = spm_slice_vol(V(3), M, V(1).dim(1:2),1);
	ty1 = A(1,1)*y1+A(1,2)*y2+A(1,3)*y3+A(1,4);
	ty2 = A(2,1)*y1+A(2,2)*y2+A(2,3)*y3+A(2,4);
	ty3 = A(3,1)*y1+A(3,2)*y2+A(3,3)*y3+A(3,4);
	VO(1) = spm_write_plane(VO(1),ty1,p);
	VO(2) = spm_write_plane(VO(2),ty2,p);
	VO(3) = spm_write_plane(VO(3),ty3,p);
    spm_progress_bar('Set',p);
end;
VO = spm_close_vol(VO);
spm_progress_bar('Clear')
return;
