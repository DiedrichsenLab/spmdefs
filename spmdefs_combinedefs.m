function [Def,mat] = spmdefs_combinedefs(Def,mat,varargin)
% function [Def,mat] = spmdefs_combinedefs(Def1,mat1,Def2,mat2,def3,mat3...)
% Combines different deformations 
% Are applied in series Def1 first, then Def1, 


Ndefs = floor(length(varargin)/2);
for i=1:Ndefs
    Def1 = Def;
    mat1 = mat;

    Def=varargin{i*2-1}; 
    mat = varargin{i*2}; 
    M    = inv(mat1);
    for j=1:size(Def{1},3)
        d0    = {double(Def{1}(:,:,j)), double(Def{2}(:,:,j)),double(Def{3}(:,:,j))};
        d{1}  = M(1,1)*d0{1}+M(1,2)*d0{2}+M(1,3)*d0{3}+M(1,4);
        d{2}  = M(2,1)*d0{1}+M(2,2)*d0{2}+M(2,3)*d0{3}+M(2,4);
        d{3}  = M(3,1)*d0{1}+M(3,2)*d0{2}+M(3,3)*d0{3}+M(3,4);
        Def{1}(:,:,j) = single(spm_sample_vol(Def1{1},d{:},[1,NaN]));
        Def{2}(:,:,j) = single(spm_sample_vol(Def1{2},d{:},[1,NaN]));
        Def{3}(:,:,j) = single(spm_sample_vol(Def1{3},d{:},[1,NaN]));
    end;
end;