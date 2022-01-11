function [Def,mat] = spmdefs_get_comp(job)
% Return the composition of a number of deformation fields.

if isempty(job),
    error('Empty list of jobs in composition');
end;
[Def,mat] = get_job(job{1});
for i=2:numel(job),
    Def1 = Def;
    mat1 = mat;
    [Def,mat] = get_job(job{i});
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
%_______________________________________________________________________

%_______________________________________________________________________
function [Def,mat] = get_job(job)
% Determine what is required, and pass the relevant bit of the
% job out to the appropriate function.

fn = fieldnames(job);
fn = fn{1};
switch fn
case {'comp'}
    [Def,mat] = spmdefs_get_comp(job.(fn));
case {'def'}
    [Def,mat] = spmdefs_get_def(job.(fn));
case {'dartel'}
    [Def,mat] = spmdefs_get_dartel(job.(fn));    
case {'sn2def'}
    [Def,mat] = spmdefs_get_sn2def(job.(fn));
case {'inv'}
    [Def,mat] = spmdefs_get_inv(job.(fn));
case {'id'}
    [Def,mat] = spmdefs_get_id(job.(fn));
otherwise
    error('Unrecognised job type');
end;
