% hidepatchborders(6);
% showpatchborders;
% setpatchborderprops('linewidth',0.5);

% extracting colorbar;
extract_colorbar = false;
if (extract_colorbar)
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
    caxis([0.000, 0.545]);
    c = colorbar(ax,'eastoutside');
    c.Label.String = '(m)';
    c.Label.VerticalAlignment = 'middle';
    % pos = get(c,'Position');
    c.Label.Position = [0.5 7];
    % ylabel(c,'(m)');
    title(c,'(m)')
    c.FontSize = 12;
    
    hec_color_range = [0.000,0.091,0.182,0.272,0.363,0.454,0.545];
    c.Ticks = hec_color_range; %create ticks
    c.TickLabels = num2cell(hec_color_range); 
    c.YAxisLocation = 'left';
    tix = c.Ticks;
    c.TickLabels = compose('%9.3f',tix);   %strictly 3dp
    ax.Visible = 'off';
    exportgraphics(gcf, 'test.eps', ...
        'ContentType', 'vector')
end

% add gauges
% add_gauges('geoclaw');

if PlotType == 4
    % hecras = load('hecras_45.m');
    hecras = load('hecras_128.dat');
    hold on
    h_hec= plot(hecras(:,1),hecras(:,2),'g',LineWidth=2);
    grid on
    % yourplots = get(gca, 'children');
    % legend(yourplots([1 2]), {'HEC-RAS', 'GeoFlood'});
    set(gca(),'Box','on');
    xlabel('Distance (m)');
    ylabel('Depth (m)') 
    xlim([0,500]);
    % title('Cross-section of depths at time T = 1 hour')
    title('')
    shg
    fig = gcf;
    fig.Units = 'inches';
    % fig.Position = [0 0 6 4]; % Adjust the width and height according to your requirements

    % load geoclaw_data
    [amrdata_ref,tref] = readamrdata_ascii(2,Frame,'./new_geoclaw_results/',false);  
    if abs(t-tref) > 1e-6
        error('times are not equal t = %g\n tref =%g',t,tref);
    end

    lstyle = {'b.','b.','b.','b.'};
    href = plotframe2ez(amrdata_ref,1,lstyle,@map1d);
    h = getlegendinfo;
    % yourplots = get(gca, 'children');
    legend([h_hec h href(1)], {'HEC-RAS', 'GeoFlood' 'GeoClaw'});
    % Save the figure as a PDF file
    filename = 'transect_90';
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
% daspect([0.50,0.502,0.505])
daspect([0.5,0.5,0.5])
% daspect([1,1,1])
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

    figsize=[4,8];
    maxlevel = 4;    %  eff. res. = 2^maxlevel
    mx = 100;
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
