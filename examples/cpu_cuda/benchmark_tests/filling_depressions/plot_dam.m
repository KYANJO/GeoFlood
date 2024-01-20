function plot_malpasset()

close all;

[ax,bx,ay,by] = plot_topo('scratch/Test6BDEM.asc','r');

% axis([957738.41,957987.1,1844520.82, 1844566.5]);
shg
axis([ax bx ay by]);
view(0.6031,13.7748)
camlight;
setviews;
view(3);

end

function [ax,bx,ay,by] = plot_topo(fname,c)
% PLOT_TOPO plots the topo from file FNAME


% c = 'w';  % Color of the topo
hold on;
[p,ax,bx,ay,by] = plot_feature(fname,c);
hold on;

%fprintf('Height at input location : %12.4f\n',hp);

daspect([1,1,.3]);
set(gca,'clipping','off')

shg

end


function [p,ax,bx,ay,by,hpout] = plot_feature(filename,cvec,xp,yp)

fid = fopen(filename, 'r');
header = textscan(fid, '%s %f', 6);
fclose(fid);

% Extract the necessary header values
ncols = header{2}(1);
nrows = header{2}(2);
xllcorner = header{2}(3);
yllcorner = header{2}(4);
dx= header{2}(5);
NODATA_value = header{2}(6);

% Read the elevation data
data = dlmread(filename, ' ', 6, 0); % Assuming values are space-separated after the header
[m,n] = size(data);
nrows = m;
ncols = n;
% Reshape the elevation data to match the grid dimensions
T = reshape(data, ncols, nrows)';
T(T== NODATA_value) = NaN; % Set NODATA values to NaN

% Create x and y coordinates for each elevation point
x = xllcorner + (0:ncols-1) * dx;
y = yllcorner + (0:nrows-1) * dx;

% --------------------------------
ax = xllcorner;
ay = yllcorner;

bx = ax + dx*(ncols-1);
by = ay + dx*(nrows-1);


[xm,ym] = meshgrid(x,y);

% nodata = find(T == nodata);
% T(nodata) = nan;
% c = T;
% c(~nodata) = 1;
% c(nodata) = 2;

p = surf(xm,ym,T);
% set(p,'cdata',c');
set(p,'edgecolor','none');

colormap([1 0 0;1 1 1]);
set(p,'cdatamapping','direct');


fprintf('Min height  : %12.4f\n',min(T(:)));
fprintf('Max height  : %12.4f\n',max(T(:)));

zlim([min(T(:)), max(T(:))]);


if (nargin > 2)
    % find height of given location (xp,yp)
    hp = interp2(xm,ym,T,xp,yp,'linear');
    if (nargout > 5)
        hpout = hp;
    end
end


end


