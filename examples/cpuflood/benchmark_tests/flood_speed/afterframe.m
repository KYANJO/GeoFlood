% hidepatchborders(6);
% showpatchborders;
% setpatchborderprops('linewidth',0.5);

colormap(flipud(parula));
% % colorbar;
tol = -0.8;
c1 = 0;
c2 = 30;
caxis([c1,c2]);

% add gauges
add_gauges('geoclaw');

fprintf('%20s %12.4e\n','qmin',qmin);
fprintf('%20s %12.4e\n','qmax',qmax);

% Axes
ax = 0;
bx = 1000;
ay = 0;
by = 2000;

axis([ax bx ay by])
daspect([0.5,0.5,0.5])
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
    figsize=[2,4];
    maxlevel = 2;    %  eff. res. = 2^maxlevel
    mx = 100;
    mi = 2;
    mj = 4;
    maxres = mi*mx*2^maxlevel;
    dpi = maxres/2
    prefix = 'flood';
    plot_tikz_fig(Frame,figsize,prefix,dpi)
    % filename = sprintf('flo_%04d.png',Frame)
    % print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
