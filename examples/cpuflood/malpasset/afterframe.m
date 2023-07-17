

colormap(flipud(parula));
% % colorbar;
tol = -0.8;
c1 = 0;
c2 = 30;
caxis([c1,c2]);

fprintf('%20s %12.4e\n','qmin',qmin);
fprintf('%20s %12.4e\n','qmax',qmax);

% Axes
ax = 953236.0;
bx = 959554.0;
ay = 1832407.25;
by = 1848572.75;

axis([ax bx ay by])
daspect([0.38 0.485 0.49])
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
% setpatchborderprops('linewidth',0.5);
    figsize=[2,4];
    maxlevel = 5;    %  eff. res. = 2^maxlevel
    mx = 32;
    mi = 1;
    mj = 2;
    maxres = mi*mx*2^maxlevel;
    dpi = maxres/figsize(1)
    prefix = 'mlp';
    plot_tikz_fig(Frame,figsize,prefix,dpi)
    % filename = sprintf('mlp_%04d.png',Frame)
    % print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
