% hidepatchborders(6);
% showpatchborders;
% setpatchborderprops('linewidth',0.5);

% colormap(flipud(parula));
% % % colorbar;
% tol = -0.8;
% c1 = 0;
% c2 = 30;
% caxis([c1,c2]);


% extracting colorbar;
% f = figure;
% ax = axes;
% mymap = [ 127 255 0 
%     0 128 0 
%     255 255 0
%     255 165 9
%     255 0 0
%     255 192 203
%     255 0 255]/255;
% Ncm = 256;
% fine_pos = linspace(0, 1, Ncm)';
% pos = linspace(0, 1, size(mymap, 1));
% mymap = interp1(pos, mymap, fine_pos);
% colormap(mymap)
% caxis([0.4, 8.0]);
% c = colorbar(ax,'eastoutside');
% ax.Visible = 'off';
% exportgraphics(gcf, 'test.eps', ...
%     'ContentType', 'vector')

% showpatchborders;
fprintf('%20s %12.4e\n','qmin',qmin);
fprintf('%20s %12.4e\n','qmax',qmax);

% % Axes
ax = -171.50000000;
bx = 1848.50000000;
ay = -56.50000000;
by = 54.50000000;

axis([ax bx ay by])
daspect([18,18,18])
% set(gca,'xtick',[]);
% set(gca,'YTick',[]);


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
    figsize=[36,2];
    maxlevel = 3;    %  eff. res. = 2^maxlevel
    mx = 54;
    mi = 18;
    mj = 1;
    maxres = mi*mx*2^maxlevel;
    dpi = maxres/figsize(1)
    prefix = 'dam_break';
    plot_tikz_fig(Frame,figsize,prefix,dpi)
    % filename = sprintf('flo_%04d.png',Frame)
    % print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;


