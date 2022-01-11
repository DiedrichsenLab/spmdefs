function spmj_defs_hotellings(kind,P1,P2,mask);
% Computes F-approximation for Hotellings T2 statistics
%   On a three-dimensional vector
%   Uses either one-sample test (against a null deformation)
%       or a twosample test

% Choices of type of Hotellings t-test

d=3;
choices={'onesample','dep','indep','indep-hetero'};
if (nargin<1)
    kind = spm_input('Choose type of test','1','b',char(choices));
end;

% Input deformations
if (nargin<2)
    P1=spm_get(inf,{'*.img','noexpand'},'Pick deformations group 1');
else
    P1=char(P1);
end;

if (nargin<3 && ~strcmp(kind,'onesample'));
    P2=spm_get(inf,{'*.img','noexpand'},'Pick deformations group 2');
else
    P2=char(P2);
end;


if (nargin<4);
    mask=spm_get([0 1],{'*.img'},'mask');
end;

if (~strcmp(kind,'onesample'));
    index = [ones(size(P1,1),1);ones(size(P2,1),1)*2];
    P=strvcat(P1,P2);
else
    P=P1;
end;

for s=1:size(P,1)
    for i=1:d
        Vi(s,i)=spm_vol(sprintf('%s,%d',deblank(P(s,:)),i));
    end;
end;

VM=spm_vol(mask);

Vo=Vi(1,1);
Vo.fname='spmF.img';
Vo=spm_create_vol(Vo);

for i=1:d
    Vd{i}= struct(		'fname',	'contr.img',...
        'dim',		Vi(1,1).dim,...
        'mat',		Vi(1,1).mat,...
        'pinfo',	Vi(1,1).pinfo,...
        'descrip',	'Mean Difference',...
        'n', i);
    Vd{i}=spm_create_vol(Vd{i});
end;

Vvar=Vi(1,1);
Vvar.fname='ResI.img';
Vvar=spm_create_vol(Vvar);

%=======================================================================
%-Computation
%=======================================================================
[n,d]   = size(Vi);			%-#images
Y   = zeros(Vo.dim(1:3));		%-result of calculations

xdim=Vi(1,1).dim(1);
ydim=Vi(1,1).dim(2);
zdim=Vi(1,1).dim(3);

%-Start progress plot
%-----------------------------------------------------------------------
spm_progress_bar('Init',Vo.dim(3),'T2','planes completed');

xords = [1:xdim]'*ones(1,ydim); xords = xords(:)';  % plane X coordinates
yords = ones(xdim,1)*[1:ydim];  yords = yords(:)';  % plane Y coordinates

%-Loop over planes computing result Y
%-----------------------------------------------------------------------
for p = 1:Vo.dim(3)
    zords   = p*ones(xdim*ydim,1)';	%-plane Z coordinates
    XM=spm_sample_vol(VM,xords,yords,zords,1);
    I=find(XM);
    Y=zeros(xdim,ydim)*NaN;
    F=[];dM=[];resI=[];
    if (length(I)>0)
    
    X=zeros(n,length(I),d);
    for j=1:d
        for i = 1:n
            data = spm_sample_vol(Vi(i,j),xords(I),yords(I),zords(I),1);
            X(i,:,j) = data';
        end
    end;

    switch (kind)
        case 'onesample'
            [F,T2,v,P]=hotT2(permute(X,[1 3 2]));
        case 'indep'
            [F,T2,v,P,dM,resI]=hotT2iho(permute(X,[1 3 2]),index);
        otherwise
            error([kind ' not implemented yet']);
    end;
    end;
    
    % Write our the images 
    Y(I) = F';
    spm_write_plane(Vo,Y,p);
    
    for (i=1:d) 
        if (~isempty(dM))
            Y(I) = dM(i,:)';    
        end;
        spm_write_plane(Vd{i},Y,p);
    end;
    
    Y(I) = resI';
    spm_write_plane(Vvar,Y,p);
    
    spm_progress_bar('Set',p);
end;


%-Write output image (uses spm_write_vol - which calls spm_write_plane)
%-----------------------------------------------------------------------
Vo = spm_close_vol(Vo);
for (i=1:d) 
    Vd{i}=spm_close_vol(Vd{i});
end;
spm_close_vol(Vvar);


sprintf('F approximation has %d, %d degrees of freedom\n',v(1),v(2));
%-End
%-----------------------------------------------------------------------
spm_progress_bar('Clear');
