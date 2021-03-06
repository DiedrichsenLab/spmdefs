function spmdefs_apply_def(Def,mat,fnames,intrp,ofnames)
% function spmdefs_apply_def(Def,mat,fnames,intrp,ofnames)
% Warp an image or series of images according to a deformation field
% INPUT: 
% Def: Deformation field data {cell}
% mat: Affine transformation matrix of the atlas (output)
% fnames: File names of the input 
% interp: interpolation 
% ofnames: outfilenames 
intrp = [intrp*[1 1 1], 0 0 0];

for i=1:size(fnames,1),
    V = spm_vol(fnames(i,:));
    M = inv(V.mat);
    [pth,nam,ext] = spm_fileparts(fnames(i,:));
    if (nargin<5)
        ofname = fullfile(pth,['w',nam,ext]);
    else 
        ofname= deblank(ofnames(i,:)); 
    end;
    Vo = struct('fname',ofname,...
                'dim',[size(Def{1},1) size(Def{1},2) size(Def{1},3)],...
                'dt',V.dt,...
                'pinfo',V.pinfo,...
                'mat',mat,...
                'n',V.n,...
                'descrip',V.descrip);
    C  = spm_bsplinc(V,intrp);
    Vo = spm_create_vol(Vo);
    for j=1:size(Def{1},3)
        d0    = {double(Def{1}(:,:,j)), double(Def{2}(:,:,j)),double(Def{3}(:,:,j))};
        d{1}  = M(1,1)*d0{1}+M(1,2)*d0{2}+M(1,3)*d0{3}+M(1,4);
        d{2}  = M(2,1)*d0{1}+M(2,2)*d0{2}+M(2,3)*d0{3}+M(2,4);
        d{3}  = M(3,1)*d0{1}+M(3,2)*d0{2}+M(3,3)*d0{3}+M(3,4);
        dat   = spm_bsplins(C,d{:},intrp);
        Vo    = spm_write_plane(Vo,dat,j);
    end;
end;
return;
