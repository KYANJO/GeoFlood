close all
clc

% Load DEM data
filename = 'scratch/test2DEM.asc'; 

% Read the header information
fid = fopen(filename, 'r');
header = textscan(fid, '%s %f', 6);
fclose(fid);

% Extract the necessary header values
ncols = header{2}(1);
nrows = header{2}(2);
xllcorner = header{2}(3);
yllcorner = header{2}(4);
cellsize = header{2}(5);
NODATA_value = header{2}(6);

% Read the elevation data
data = dlmread(filename, ' ', 6, 0); % Assuming values are space-separated after the header

% Reshape the elevation data to match the grid dimensions
elevation = reshape(data, ncols, nrows)';
elevation(elevation == NODATA_value) = NaN; % Set NODATA values to NaN

% Create x and y coordinates for each elevation point
x = xllcorner + (0:ncols-1) * cellsize;
y = yllcorner + (0:nrows-1) * cellsize;

% axis([0 2000 0 2000])

% Create a grid for the elevation data
[X, Y] = meshgrid(x, y);

% Plot the DEM in 3D
figure;

surf(X, Y, elevation, 'FaceColor', 'interp', 'EdgeColor', 'none');
hold on;

% Add contour lines
contour3(X, Y, elevation, 20, 'k'); % Change the number 20 to adjust the number of contour lines

colormap(white); % Use the parula colormap for the surface
% colorbar;
% 
% draw the square boundary
hold on
s = linspace(0,2000,100);
plot3(0*s,s,interp2(X,Y,elevation,0*s,s,'linear'),'m','LineWidth',2);
plot3(s,0*s,interp2(X,Y,elevation,s,0*s,'linear'),'m','LineWidth',2);
plot3(0*s+2000,s,interp2(X,Y,elevation,0*s+2000,s,'linear'),'m','LineWidth',2);
plot3(s,0*s+2000,interp2(X,Y,elevation,s,0*s+2000,'linear'),'m','LineWidth',2);


% plot3([0,0],[0,2000],interp2(X, Y, elevation, [0,0],[0,2000], 'linear'),'m','LineWidth',1);
% plot3([0,2000],[2000,2000],interp2(X, Y, elevation, [0,2000],[2000,2000], 'linear'),'m','LineWidth',1);
% plot3([2000,2000],[2000,0],interp2(X, Y, elevation, [2000,2000],[2000,0], 'linear'),'m','LineWidth',1);
% plot3([2000,0],[0,0],interp2(X, Y, elevation, [2000,0],[0,0], 'linear'),'m','LineWidth',1);

% Define the boundary conditions
startRow = 1900;
startCol = 0;
endRow = 2000;
endCol = 0;

% Plot the line on the DEM
line_x = [startCol, endCol];
line_y = [startRow, endRow];
line_z = interp2(X, Y, elevation, line_x, line_y, 'linear');
hold on;
plot3(line_x, line_y, line_z, 'b', 'LineWidth', 8);

% drawing control points 1,4,10,12
hold on;
plot_points(X,Y,250,250, 'Point 1',elevation);
plot_points(X,Y,250, 1750, 'Point 4',elevation);
plot_points(X,Y,1250,750, 'Point 10',elevation);
plot_points(X,Y,1250,1750, 'Point 12',elevation);

% Adjust the viewing angle and perspective
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Elevation (m)');
% title('3D Plot of DEM');

% camlight;
zlim([7,11]);
axis on
setviews;
% view(3);
% [az,el] = view;
% view([-29.1972,68.2975]);
light               % create a light
lighting gouraud    % preferred method for lighting curved surfaces
shg
% axis equal
% daspect([1 1 0.0015]); view([-29.1972,68.2975]);
daspect([1 1 1])

% rotate3d on

% axis([x(1) x(end) y(1) y(end) min(elevation(:)) max(elevation(:))]);
set(gca, 'FontWeight', 'bold');
% saveas(gcf, 'depressions.pdf', 'pdf');
% Adjust the size and position of the figure
fig = gcf;
fig.Units = 'inches';
% fig.Position = [0 0 6 4]; % Adjust the width and height according to your requirements

% Save the figure as a PDF file
filename = 'figure_flat';
exportgraphics(fig, [filename '.jpeg'], 'ContentType', 'vector', 'Resolution', '700');

function [a] = plot_points(X,Y,x,y,point, elevation)

%     p1 = elevation(y,x); %(y,x)
    p1 = interp2(X, Y, elevation, x, y);
    scatter3(x,y,p1+10,'red','filled')
    text(x, y, p1 + 75, point, 'FontSize', 12, 'FontWeight', 'bold');

end