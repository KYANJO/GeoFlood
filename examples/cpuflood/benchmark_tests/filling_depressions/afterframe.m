% hidepatchborders(6);
% showpatchborders;
% setpatchborderprops('linewidth',0.5);

% colormap(flipud(parula));
% % % colorbar;
% tol = -0.8;
% c1 = 0;
% c2 = 30;
% caxis([c1,c2]);
% 
% % add gauges
% add_gauges('geoclaw');

% extracting colorbar;
f = figure;
ax = axes;
mymap = [ 127 255 0 
    0 128 0 
    255 255 0
    255 165 9
    255 0 0
    255 192 203
    255 0 255]/255;
Ncm = 256;
fine_pos = linspace(0, 1, Ncm)';
pos = linspace(0, 1, size(mymap, 1));
mymap = interp1(pos, mymap, fine_pos);
colormap(mymap)
caxis([7.95, 10.36]);
c = colorbar(ax,'eastoutside');
hec_color_range = [7.95,8.35,8.76,9.16,9.56,9.96,10.36];
c.Ticks = hec_color_range; %create ticks
c.TickLabels = num2cell(hec_color_range); 
ax.Visible = 'off';
exportgraphics(gcf, 'test.eps', ...
    'ContentType', 'vector')

% colorbar 
% showpatchborders;
fprintf('%20s %12.4e\n','qmin',qmin);
fprintf('%20s %12.4e\n','qmax',qmax);

% Axes
ax = 0;
bx = 2000;
ay = 0;
by = 2000;

axis([ax bx ay by])
daspect([1,1,1])
set(gca,'xtick',[]);
set(gca,'YTick',[]);


% title(sprintf('malpasset Dam (%d) : t = %.2f (%.2f,%.2f)',Frame,t,qmin,qmax),'fontsize',18);
title('')
ylabel('')
xlabel('')

NoQuery = 0;
prt = true;
MaxFrames = 1000;
if (prt)
    % hidepatchborders(6);
    % setpatchborderprops('linewidth',0.5);
    % % showpatchborders(1:5)
    figsize=[4,4];
    maxlevel = 2;    %  eff. res. = 2^maxlevel
    mx = 50;
    mi = 4;
    mj = 8;
    maxres = mi*mx*2^maxlevel;
    dpi = maxres/figsize(1)
    prefix = 'filling';
    plot_tikz_fig(Frame,figsize,prefix,dpi)
    % filename = sprintf('flo_%04d.png',Frame)
    % print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;


