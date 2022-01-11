function def_plot_vecfield(fname,varargin);
% Plots a deformation vector field as 2-d vector field 
% function def_plot_vecfield(name,varargin);
% varargin: 
% 'orientation' (default coronal)
% 'slice' (default 40)
% 'spacing' (default 5)
% 'scale' (default 0 for automatic) 
% Flags: 
% 'Ticks'
orientation='y';
underlay=[];
mask=[];
affine=eye(4); 
slice=40;
spacing=5;
Ticks=1;
scale=0; % 
vararginoptions(varargin,{'','','orientation','slice','spacing','scale','underlay','mask','affine'},{'Ticks'});

% Get the vector field 
if (nargin<1 | isempty(fname))
    fname=spm_select(1,'image','Get vector field ');
end;

% Read in the vector field 
for i=1:3
    V(i)=spm_vol(sprintf('%s,1,%d',fname,i));
end;
% Get the corner of the image as a bounding box 
[corx,cory,corz]=spmj_affine_transform([1 V(1).dim(1)],[1 V(1).dim(2)],[1 V(1).dim(3)],V(1).mat); 

% Generate Meshgrid over the volume in world space 
% Be careful as Meshgrid puts the first dimension into the columns 
% and second dimnesion into rows.... Image is exactly opposite 
switch(orientation)
    case 'z'
        [I{2},I{1},I{3}]=meshgrid([min(cory):spacing:max(cory)],[min(corx):spacing:max(corx)],slice);
        [Ia{2},Ia{1},Ia{3}]=meshgrid([min(cory):1:max(cory)],[min(corx):1:max(corx)],slice);
        indx=[1 2];
    case 'y'
        [I{2},I{1},I{3}]=meshgrid(slice,[min(corx):spacing:max(corx)],[min(corz):spacing:max(corz)]);
        [Ia{2},Ia{1},Ia{3}]=meshgrid(slice,[min(corx):1:max(corx)],[min(corz):1:max(corz)]);
        indx=[1 3];
    case 'x' 
        [I{2},I{1},I{3}]=meshgrid([min(cory):spacing:max(cory)],slice,[min(corz):spacing:max(corz)]);
        [Ia{2},Ia{1},Ia{3}]=meshgrid([min(cory):1:max(cory)],slice,[min(corz):1:max(corz)]);
        indx=[2 3];
    otherwise 
        error('Orientation has to be x,y,z');
end;

% Squeeze
for i=1:3
    I{i}=squeeze(I{i});
    Ia{i}=squeeze(Ia{i});
end;

% Sample the underlay from the anatomical image 
if (~isempty(underlay))
    VA=spm_vol(underlay);
    [Iav{1},Iav{2},Iav{3}]=spmj_affine_transform(Ia{1},Ia{2},Ia{3},inv(VA.mat)); % Transfer into voxel space of anatomy
    A=spm_sample_vol(VA,Iav{1},Iav{2},Iav{3},1);          
    sc(1)=min(A(:));
    sc(2)=max(A(:));
    image(Ia{indx(1)}(:,1),Ia{indx(2)}(1,:),(A'-sc(1))./(sc(2)-sc(1))*255);
    C=[0:1/254:1]';C=[C C C];
    colormap(C);
       set(gca,'YDir','normal');
    hold on;
end;

% Sample Vector map
[Iv{1},Iv{2},Iv{3}]=spmj_affine_transform(I{1},I{2},I{3},inv(V(1).mat)); % Transfer into voxel space of vector 
for i=1:3 
    Y{i}=spm_sample_vol(V(i),Iv{1},Iv{2},Iv{3},1);
end;

if (~isempty(mask))
    VM=spm_vol(mask);
    [Iv{1},Iv{2},Iv{3}]=spmj_affine_transform(I{1},I{2},I{3},inv(VM.mat)); % Transfer into voxel space of mask
    XM=spm_sample_vol(VM,I{1},I{2},I{3},1);
    Y{1}(XM==0)=NaN;
    Y{2}(XM==0)=NaN;
    Y{3}(XM==0)=NaN;
end;


% Plot the vector map
if (scale==0)
    quiver(I{indx(1)},I{indx(2)},Y{indx(1)},Y{indx(2)});
else 
    quiver(I{indx(1)},I{indx(2)},Y{indx(1)}*scale,Y{indx(2)}*scale,0);
end;

if (~Ticks);
    set(gca,'XTick',[]);
    set(gca,'YTick',[])
end;
set(gca,'Box','off');axis equal;
hold off;