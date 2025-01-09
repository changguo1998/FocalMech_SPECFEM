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
[Yq1, Xq1, Zq1] = meshgrid([125], xs, zs);
[Yq2, Xq2, Zq2] = meshgrid(ys, [184.5], zs);
[Yq3, Xq3, Zq3] = meshgrid(ys, xs, [1.0]);

vari = vp;

cmap = [1 1 1; interp1([0, 50, 100], [240 0 86; 250 250 250; 40 40 238]./255, 0:100)];
%%
figure('Position', [100, 100, 1200, 400]);
S1 = interp3(Y, X, Z, vari, Yq1, Xq1, Zq1);
ax1 = axes;
imagesc([xs(1), xs(end)], [zs(1), zs(end)], reshape(S1, [size(S1, 1), size(S1, 3)])'); hold on;
cb1 = colorbar;
minvari = min(S1(S1>0.4))-0.1;
maxvari = max(S1(:))+0.1;
if minvari == maxvari
    minvari = maxvari - 0.1;
    maxvari = maxvari + 0.1;
end
caxis([5, maxvari]);
colormap(cmap);
plot(184.5, -2.19, 'k^', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
xlabel('北 (km)');
ylabel('深度 (km)');
ax1.XLabel.FontSize = 16;
ax1.YLabel.FontSize = 16;
ax1.FontSize = 14;
cb1.Label.String = '速度 (km/s)';
cb1.Label.FontSize = 16;

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('-painters', '-dpdf', '-r600', '/data/chuan-dian-event_lt4.0/src/paper_plot/YNYOS_FDmodel_sliceNS', '-bestfit');
exportgraphics(gcf, '/data/chuan-dian-event_lt4.0/src/paper_plot/YNYOS_FDmodel_sliceNS.png', 'Resolution', 300);
%%
figure('Position', [100, 100, 1200, 400]);
S2 = interp3(Y, X, Z, vari, Yq2, Xq2, Zq2);
ax2 = axes;
imagesc([ys(1), ys(end)], [zs(1), zs(end)], reshape(S2, [size(S2, 2), size(S2, 3)])');hold on;
cb2 = colorbar;
minvari = min(S2(S2>0.4))-0.1;
maxvari = max(S2(:))+0.1;
if minvari == maxvari
    minvari = maxvari - 0.1;
    maxvari = maxvari + 0.1;
end
caxis([5, maxvari]);
colormap(cmap);
plot(125, -2.19, 'k^', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
xlabel('东 (km)');
ylabel('深度 (km)');
ax2.XLabel.FontSize = 16;
ax2.YLabel.FontSize = 16;
ax2.FontSize = 14;
cb2.Label.String = '速度 (km/s)';
cb2.Label.FontSize = 16;



fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('-painters', '-dpdf', '-r600', '/data/chuan-dian-event_lt4.0/src/paper_plot/YNYOS_FDmodel_sliceWE', '-bestfit');
exportgraphics(gcf, '/data/chuan-dian-event_lt4.0/src/paper_plot/YNYOS_FDmodel_sliceWE.png', 'Resolution', 300);

%%

figure('Position', [100, 100, 800, 700]);
S3 = interp3(Y, X, Z, vari, Yq3, Xq3, Zq3);
ax3 = axes;
imagesc([ys(1), ys(end)], [xs(1), xs(end)], S3); hold on;
cb3 = colorbar;
minvari = min(S3(S3>0.4));
maxvari = max(S3(:));
if minvari == maxvari
    minvari = maxvari - 0.1;
    maxvari = maxvari + 0.1;
end
caxis([minvari, maxvari]);
colormap(ax3, cmap(2:end, :));
plot([0, 400], [184.5, 184.5], '-', 'LineWidth', 2, 'Color', [1, 1, 1]./2);
plot([125, 125], [0, 400], '-', 'LineWidth', 2, 'Color', [1, 1, 1]./2);
plot(125, 184.5, 'k^', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
set(ax3, 'YDir', 'normal');
xlabel('东 (km)');
ylabel('北 (km)');
ax3.XLabel.FontSize = 16;
ax3.YLabel.FontSize = 16;
ax3.FontSize = 14;
cb3.Label.String = '速度 (km/s)';
cb3.Label.FontSize = 16;


fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('-painters', '-dpdf', '-r600', '/data/chuan-dian-event_lt4.0/src/paper_plot/YNYOS_FDmodel_sliceDep1km', '-bestfit');
exportgraphics(gcf, '/data/chuan-dian-event_lt4.0/src/paper_plot/YNYOS_FDmodel_sliceDep1km.png', 'Resolution', 300);

% Ss = slice(Y, X, Z, vari, [220], [380], [1], 'nearest');
% for i = 1:length(Ss)
%     Ss(i).EdgeColor = 'none';
% end
% colorbar;
% xlabel('East');
% ylabel('North');
% zlabel('Depth');
% colormap(interp1([1; 33; 65], [1.0, 0.0, 0.0; 1.0, 1.0, 1.0; 0.0, 0.0, 1.0].*0.9, 1:65));
% axis equal;
% minvari = min(vari(vari>0.4));
% maxvari = max(vari(:));
% if minvari == maxvari
%     minvari = maxvari - 0.1;
%     maxvari = maxvari + 0.1;
% end
% caxis([minvari, maxvari]);
% set(gca, 'ZDir', 'reverse');