clc;
clear all;
close all;

dir = 'G:\cell cycle and differentiation of stem cell\mouse retinal development\landscape\cosine\0.1_0.1\';

dir = 'G:\HSC\HSC_landscape\multiprocess\';
% Load the landscape data
Xgrid = csvread([dir,'Xgrid.csv']);
Ygrid = csvread([dir,'Ygrid.csv']);
PotU = csvread([dir,'pot_U.csv']);

% Create a meshgrid for plotting
% [X, Y] = meshgrid(Xgrid, Ygrid);
U = PotU';
% Plot the landscape
figure;
surf(Xgrid, Ygrid, U);
colormap('turbo');
colormap('jet');
alpha(1);
shading interp;
lighting gouraud;
xlabel('UMAP1');
ylabel('UMAP2');
zlabel('Potential');
set(gca, 'FontName', 'Arial');
set(gca,'FontSize',16, 'LabelFontSizeMultiplier', 1, 'TitleFontSizeMultiplier', 1);
pbaspect([1 0.9 0.8]);
set(gca,'TickDir', 'out', 'TickLength', [0.02 0.02])
set(gca, 'LineWidth', 2, 'Color', [0 0 0])
set(gca, 'XColor', [0.00 0.00 0.00])
set(gca, 'YColor', [0.00 0.00 0.00])
set(gca, 'ZColor', [0.00 0.00 0.00])
xlim([-10, 16])  %mouse
ylim([-8, 15])  %mouse
zlim([4, 20])  %mouse
xticks([-10, 0, 10]);   
yticks([-5, 5, 15]);
%yticks([-7, 0, 7, 14]) %mouse
zticks([10, 20]) %mouse

xlim([ 7.8, 20]) %hsc
ylim([ 0, 17.01153336]) %hsc
zlim([ 5, 16.3]) %hsc
xticks([8, 12, 16, 20]) %HSC
yticks([0, 5, 10, 15]) %HSC
zticks([5, 10, 15]) %HSC
set(gca, 'color', 'white')
%colorbar('Ticks', [6, 10, 14, 18], 'TickLength', 0.02, 'TickDirection', 'out', 'FontSize', 20, 'Color', [0 0 0], 'LineWidth', 2)
colorbar('Ticks', [8, 10, 12, 14, 16], 'TickLength', 0.02, 'TickDirection', 'out', 'FontSize', 20, 'Color', [0 0 0], 'LineWidth', 2)
box on 

hold on
mesh(Xgrid(1:6:end, 1:6:end), Ygrid(1:6:end, 1:6:end), U(1:6:end, 1:6:end)+0.2, 'LineStyle', '-', 'LineWidth', 0.7, 'EdgeColor', 'k', 'FaceColor', 'none')
hold on;


% Load and plot each path
paths = {'Progenitor_PR.csv', 'Progenitor_AC-HC.csv', 'Progenitor_RGC.csv', ...
         'PR_Progenitor.csv', 'PR_AC-HC.csv', 'PR_RGC.csv', ... 
         'AC-HC_Progenitor.csv', 'AC-HC_PR.csv', 'AC-HC_RGC.csv', ... 
         'RGC_Progenitor.csv', 'RGC_PR.csv', 'RGC_AC-HC.csv'};

paths = {'Bas_Ery.csv', 'Bas_HSC.csv', 'Bas_Meg.csv', 'Bas_Mon.csv', 'Bas_Neu.csv', ... % and other paths
         'Ery_Bas.csv', 'Ery_HSC.csv', 'Ery_Meg.csv', 'Ery_Mon.csv', 'Ery_Neu.csv', ...
         'HSC_Bas.csv', 'HSC_Ery.csv', 'HSC_Meg.csv', 'HSC_Mon.csv', 'HSC_Neu.csv', ...
         'Meg_Bas.csv', 'Meg_Ery.csv', 'Meg_HSC.csv', 'Meg_Mon.csv', 'Meg_Neu.csv', ...
         'Mon_Bas.csv', 'Mon_Ery.csv', 'Mon_HSC.csv', 'Mon_Meg.csv', 'Mon_Neu.csv', ...
         'Neu_Bas.csv', 'Neu_Ery.csv', 'Neu_HSC.csv', 'Neu_Mon.csv', 'Neu_Meg.csv'};

dir1 = 'G:\HSC\HSC_landscape\D_xyGridSpacing_numTimeSteps_dt_Tragrid\0.2_1_2000000_0.2_200\minmax_fantan\';

for i = 1:length(paths)
    path = csvread([dir1,paths{1,i}],1,0);
    U_values = interpolateU(Xgrid, Ygrid, U, path);
    plot3(path(:,1), path(:,2), U_values+0.25, 'LineWidth', 2.5);
end

hold on;
dir1 = 'G:\cell cycle and differentiation of stem cell\mouse retinal development\figure\plot_gradient_path\';
gradient_paths = {'path_Pro_PR.mat', 'path_Pro_AC.mat', 'path_Pro_RGC.mat', ...
         'path_PR_Pro.mat', 'path_PR_AC.mat', 'path_PR_RGC.mat', ... 
         'path_AC_Pro.mat', 'path_AC_PR.mat', 'path_AC_RGC.mat', ... 
         'path_RGC_Pro.mat', 'path_RGC_PR.mat', 'path_RGC_AC.mat'};
for i = 1:length(gradient_paths)
    path_data = load([dir1,gradient_paths{1,i}]);
    path = path_data.path;
    path = path';
    U_values = interpolateU(Xgrid, Ygrid, U, path);
    plot3(path(:,1), path(:,2), U_values+0.25, 'w', 'LineWidth', 2.5);
end
hold on;

saveas( figure(1), [dir, 'landscape_path3.fig']); 
print(figure(1), '-r600', '-dpdf', [dir, 'landscape_path3.pdf']);
