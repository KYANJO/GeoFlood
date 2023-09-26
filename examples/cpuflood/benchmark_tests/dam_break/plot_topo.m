close all
clc

% Load DEM data
filename = 'scratch/Test6BDEM.asc'; 

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
[m,n] = size(data);
nrows = m;
ncols = n;
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
colormap(parula); % You can change the colormap to your preference
colorbar;
xlabel('X');
ylabel('Y');
zlabel('Elevation');
title('3D Plot of DEM');

% Adjust the viewing angle and perspective
view(3);

