% hidepatchborders(6);
% showpatchborders(1:5);
% setpatchborderprops('linewidth',0.5);

colormap(flipud(parula));
% % colorbar;
tol = -0.8;
c1 = 0;
c2 = 30;
caxis([c1,c2]);

% add gauges
add_gauges('forestclaw');

fprintf('%20s %12.4e\n','qmin',qmin);
fprintf('%20s %12.4e\n','qmax',qmax);

% Axes
ax = 0;
bx = 2000;
ay = 0;
by = 2000;

axis([ax bx ay by])
daspect([0.4 0.4 0.4])
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
%     hidepatchborders(6);
%     setpatchborderprops('linewidth',0.5);
    showpatchborders;
    figsize=[4,4];
    maxlevel = 1;    %  eff. res. = 2^maxlevel
    mx = 40;
    mi = 5;
    mj = 5;
    maxres = mi*mx*2^maxlevel;
    dpi = maxres/4
    prefix = 'plot';
    plot_tikz_fig(Frame,figsize,prefix,dpi)
    filename = sprintf('plot_%04d.png',Frame)
    print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
