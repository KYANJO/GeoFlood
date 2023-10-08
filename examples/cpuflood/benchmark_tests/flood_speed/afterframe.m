% hidepatchborders(6);
% showpatchborders;
% setpatchborderprops('linewidth',0.5);

% colormap(flipud(parula));
% mymap = [ 0.498 1 0 
%     0 0.5 0 
%     1 1 0 
%     1 0.647 0 
%     1 0 0 
%     1 0.752 0.796 
%     1 0 1];
% 
% pos = linspace(0, 1, size(mymap, 1));
% 
% % Create a fine-grained colormap using interpolation
% fine_pos = linspace(0, 1, 256);
% fine_mymap = interp1(pos, mymap, fine_pos);
% % Set the stretched colormap
% colormap(fine_mymap);
% colorbar;

% tol = -0.8;
% c1 = 0;
% c2 = 30;
% caxis([c1,c2]);

% add gauges
% add_gauges('geoclaw');

if PlotType == 4
    hecras = load('hecras_45.m');
    hold on
    plot(hecras(:,1),hecras(:,2),'r',LineWidth=2);
    grid on
    yourplots = get(gca, 'children');
    legend(yourplots([1 2]), {'HEC-RAS', 'GeoFlood'});
    xlabel('Distance (m)');
    ylabel('Depth (m)') 
    xlim([0,500]);
    title('Cross-section of depths at time T = 1 hour')
    shg
    fig = gcf;
    fig.Units = 'inches';
    % fig.Position = [0 0 6 4]; % Adjust the width and height according to your requirements

    % Save the figure as a PDF file
    filename = 'transect_1_5';
    exportgraphics(fig, [filename '.jpeg'], 'ContentType', 'vector', 'Resolution', '700');
    return 

end


fprintf('%20s %12.4e\n','qmin',qmin);
fprintf('%20s %12.4e\n','qmax',qmax);

% Axes
ax = 0;
bx = 1000;
ay = 0;
by = 2000;

axis([ax bx ay by])
daspect([0.50,0.502,0.505])
% daspect([1,1,1])
% set(gca,'xtick',[]);
% set(gca,'YTick',[]);



% title(sprintf('malpasset Dam (%d) : t = %.2f (%.2f,%.2f)',Frame,t,qmin,qmax),'fontsize',18);
title('')
ylabel('')
xlabel('')


NoQuery = 0;
prt = false;
MaxFrames = 1000;
if (prt)
    % hidepatchborders(6);
    % setpatchborderprops('linewidth',0.5);
    % % showpatchborders(1:5)

    figsize=[4,8];
    maxlevel = 4;    %  eff. res. = 2^maxlevel
    mx = 50;
    mi = 2;
    mj = 4;
    maxres = mi*mx*2^maxlevel;
    dpi = maxres/figsize(1)
    prefix = 'flood';
    plot_tikz_fig(Frame,figsize,prefix,dpi)
    % colormap("parula");
    % colorbar;
    filename = sprintf('flo_%04d.png',Frame)
    print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
