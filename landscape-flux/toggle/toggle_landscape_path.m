clc;
clear all;
close all;

dir = 'G:\HSC\HSC_landscape\toggle\';

% Load the landscape data
traj = textread([dir,'landscape.txt']);
P = traj(:,3);
P = P/sum(P);
U= -log(P);
U(~isfinite(U))=18;
U(U>15)=15;

x=0:6/199:6;
y=0:6/199:6;

[xx,yy]=meshgrid(x,y);
UU=reshape(U,[200,200]);

surf(xx, yy, UU);
colormap('turbo');
alpha(1);
shading interp;
lighting gouraud;
xlabel('X');
ylabel('Y');
zlabel('Potential');
set(gca, 'FontName', 'Arial');
set(gca,'FontSize',20, 'LabelFontSizeMultiplier', 1, 'TitleFontSizeMultiplier', 1);
pbaspect([1 0.9 0.8]);
set(gca,'TickDir', 'out', 'TickLength', [0.02 0.02])
set(gca, 'LineWidth', 2, 'Color', [0 0 0])
set(gca, 'XColor', [0.00 0.00 0.00])
set(gca, 'YColor', [0.00 0.00 0.00])
set(gca, 'ZColor', [0.00 0.00 0.00])
set(gca, 'color', 'white')
box on 
hold on

% Load and plot each path
paths = {'path_coordinates_1.csv', 'path_coordinates_2.csv'};

for i = 1:length(paths)
    path = csvread([dir,paths{1,i}],1,0);
    U_values = interpolateU(xx, yy, UU, path);
    plot3(path(:,1), path(:,2), U_values+0.25, 'LineWidth', 2.5);
end
hold on;

% Define the starting and destination stable states
A_state = [4.73087644; 0.20759548]; 
B_state = [0.22273542; 4.73670019];

%%%%%% gradient path  %%%%%%%%%
format long;
global Eeff;

X = xx;
Y = yy;
U = UU;

parameters.Dim    = 2;   
parameters.N      = 50;
parameters.lambda = 100;

parameters.Xi = [4.73087644; 0.20759548]; %A
parameters.Xf = [0.22273542; 4.73670019]; %B

Eeff = 0;

x0 = zeros(parameters.Dim,parameters.N);
for j = 1:parameters.N
        x0(:,j) = parameters.Xi+ j*(parameters.Xf - parameters.Xi)/(parameters.N+1);    
end


shj0 = SHJ_WKB(X,Y,U,x0,parameters);  

options  =  optimset( 'Display',      'iter',                                 ...
                      'MaxFunEvals',  1000 * parameters.Dim * parameters.N,  ...
                      'MaxIter',      300 ,                                  ...%);
                      'PlotFcns',     @TPlot           );
%                   
% options  =  saoptimset('simulannealbnd')
% options = saoptimset('PlotFcns',{@saplotbestx,...
%                 @saplotbestf,@saplotx,@saplotf});              
                              
tic
[xxx, fval, exitflag]  =  fminunc(@(x) SHJ_WKB(X,Y,U,x,parameters), x0, options);
% [xxx, fval, exitflag]  =  simulannealbnd(@(x) SHJ_WKB(x,parameters), x0,[],[],options);
toc

fval 
exitflag 

path = [double(parameters.Xi) xxx double(parameters.Xf)];
file=sprintf('path_A_B.mat');
file=sprintf('path_B_A.mat');




saveas( figure(1), [dir, 'figure/landscape_path.fig']); 
print(figure(1), '-r600', '-dpdf', [dir, 'figure/landscape_path.pdf']);
