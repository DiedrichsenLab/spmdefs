function spmdefs_apply_def(Def,mat,fname,intrp,ofname)
% function spmdefs_apply_def(Def,mat,fnames,intrp,ofnames)
% Warp an image or series of images according to a deformation field
% INPUT: 
% Def: Deformation field data {cell}
% mat: Affine transformation matrix of the atlas (output)
% fname: File name of the input 
% interp: interpolation 
% ofname: outfilename 

intrp = [intrp*[1 1 1], 0 0 0];

[pth,nam,ext] = spm_fileparts(fname);
if (nargin<5)
    ofname = fullfile(pth,['w',nam,ext]);
else 
    ofname= deblank(ofname); 
end;

V = spm_vol(fname);
    
for i=1:length(V)
    M = inv(V(i).mat);
    Vo = struct('fname',ofname,...
                'dim',[size(Def{1},1) size(Def{1},2) size(Def{1},3)],...
                'dt',V(i).dt,...
                'pinfo',V(i).pinfo,...
                'mat',mat,...
                'n',V(i).n,...
                'descrip',V(i).descrip);
    C  = spm_bsplinc(V(i),intrp);
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
