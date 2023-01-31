function plot_malpasset()

close all;

[ax,bx,ay,by] = plot_topo('scratch/Malpasset/malpasset_domaingrid_20m_nolc.topotype2','g');
% [ax,bx,ay,by] = plot_topo('scratch/Malpasset/malpasset_resevoir_5m_nolc.topotype2','b');
% [ax,bx,ay,by] = plot_topo('scratch/Malpasset/malpasset_grid4_2m_nolc.topotype2','b'); 
% [ax,bx,ay,by] = plot_topo('scratch/Malpasset/malpasset_grid3_2m_nolc.topotype2','b');
% [ax,bx,ay,by] = plot_topo('scratch/Malpasset/malpasset_grid2_1m_nolc.topotype2','r');
% [ax,bx,ay,by] = plot_topo('scratch/Malpasset/malpasset_damapproach_1m_nolc.topotype2','w');

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


function [p,ax,bx,ay,by,hpout] = plot_feature(fname,cvec,xp,yp)

fid = fopen(fname);

ncols = fscanf(fid,'%d',1); fscanf(fid,'%s',1);
nrows = fscanf(fid,'%d',1); fscanf(fid,'%s',1);
xll = fscanf(fid,'%g',1);   fscanf(fid,'%s',1);
yll = fscanf(fid,'%g',1);   fscanf(fid,'%s',1);
dx = fscanf(fid,'%g',1);    fscanf(fid,'%s',1);
nodata = fscanf(fid,'%g',1); fscanf(fid,'%s',1);
T = fscanf(fid,'%g',nrows*ncols);
fclose(fid);

% --------------------------------
ax = xll;
ay = yll;

bx = ax + dx*(ncols-1);
by = ay + dx*(nrows-1);

x = linspace(ax,bx,ncols);
y = linspace(ay,by,nrows);

T = reshape(T,ncols,nrows);
T = fliplr(T);
% T = reshape(T,nrows,ncols)';
% T = fliplr(T);


[xm,ym] = meshgrid(x,y);

nodata = find(T == nodata);
T(nodata) = nan;
c = T;
c(~nodata) = 1;
c(nodata) = 2;

p = surf(xm,ym,T');
set(p,'cdata',c');
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


