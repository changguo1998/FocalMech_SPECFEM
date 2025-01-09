clc;close all;
% modelfilename='../test/model_homo.bin';
% modelfilename='/data/glibshelltest/sem/shotN/DATA/meshfem3D_files/fdmodel.bin';
% modelfilename='/data/glibshelltest/model_3layer.bin';
% modelfilename = '/data/chuandian_sim/chuandian.model';
% modelfilename = '/home/guochang/Projects/chuandian_simulate/chuandian.model';
% modelfilename = '/data/chuandian_sim/YN.YAJ.model';
% modelfilename = '/data/testsemshell/synevent/model/changning.bin';
if ~exist('modelfilename', 'var')
    [modelf, modelfilepath] = uigetfile('*.bin');
    modelfilename = fullfile(modelfilepath, modelf);
end

fid = fopen(modelfilename, 'r');
nf = fread(fid, 1, 'uint32');
v0 = fread(fid, nf, 'float32');
dv = fread(fid, nf, 'float32');
np = fread(fid, nf, 'int32');
nm = fread(fid, 1, 'int32');
table = reshape(fread(fid, nf*nm, 'float32'), [nf, nm]);
nx = fread(fid, 1, 'int32');
ny = fread(fid, 1, 'int32');
nz = fread(fid, 1, 'int32');
dx = fread(fid, 1, 'float32');
dy = fread(fid, 1, 'float32');
dz = fread(fid, 1, 'float32');
zt = fread(fid, 1, 'float32');
index = reshape(fread(fid, nx*ny*nz, 'int32'), [nz, ny, nx]);
fclose(fid);
vp = zeros(nx,ny,nz);
vs = zeros(nx,ny,nz);
rho = zeros(nx,ny,nz);
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            idx = index(k,j,i);
            if idx > 0
                vp(i,j,k) = table(1,idx);
                vs(i,j,k) = table(2,idx);
                rho(i,j,k) = table(3,idx);
            else
                vp(i,j,k) = 0.34;
                vs(i,j,k) = 0.001;
                rho(i,j,k) = 0.00129;
            end
        end
    end
end
xs = (0:nx-1)*dx;
ys = (0:ny-1)*dy;
zs = (0:nz-1)*dz+zt;


[Y, X, Z] = meshgrid(ys, xs, zs);

vari = vp;

Ss = slice(Y, X, Z, vari, [125], [184.5], [1], 'nearest');
for i = 1:length(Ss)
    Ss(i).EdgeColor = 'none';
end
colorbar;
xlabel('East');
ylabel('North');
zlabel('Depth');
% colormap(interp1([1; 33; 65], [1.0, 0.0, 0.0; 1.0, 1.0, 1.0; 0.0, 0.0, 1.0].*0.9, 1:65));
% colormap([1.0, 1.0, 1.0; parula]);
axis equal;
minvari = min(vari(vari>0.4));
maxvari = max(vari(:));
if minvari == maxvari
    minvari = maxvari - 0.1;
    maxvari = maxvari + 0.1;
end
caxis([minvari, maxvari]);
set(gca, 'ZDir', 'reverse');