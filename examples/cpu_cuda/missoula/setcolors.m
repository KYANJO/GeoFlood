function setcolors(p,x,y,z,q)
    %
    % SETCOLORS specifies user-deqd mapping between q data and color map
    %
    %      User-deqd routine which maps q values to the colormap.  To
    %      indicate that this routine should be called, the user should set
    %      the flag 'UserColorMapping' to 1 in SETPLOT2 or SETPLOT3 files.
    %
    %      The syntax for this routine is
    %
    %              setcolors(p,x,y,z,q);
    %
    %      where p is the handle to the patch whose colors are being set,
    %      (x,y,z) are the Cartesian locations of cell centers with corresponding
    %      data values in q.  If a mapped grid or manifold is being plotted,
    %      (x,y,z) are given in the Cartesian locations, NOT the physical
    %      locations of cell centers.
    %
    %      Possible uses for this routine include :
    %
    %             -- Use different colormaps in different domain regions
    %             -- Mask out embedded boundary regions.
    %             -- Visualize parallel partitions
    %
    %      To visualize under and over shoots, see UNDEROVER.
    %
    %      Example : Use the 'spring' colormap for the left half of a region
    %      and the 'winter' color map for the right half.
    %
    %              c1 = spring;
    %              c2 = winter;
    %
    %              cm_new = [c1; c2];
    %
    %              n1 = length(c1);
    %              n2 = length(c2);
    %
    %              % Map q(m1) to indices [1,n1] and q(m2) to indices [n1+1,n1+n2]
    %              idx = 0*q + nan;   % To make sure we don't leave any q values out.
    %
    %              m1 = x > 0.5;      % Left and right halves of computational domains
    %              m2 = x <= 0.5;
    %
    %              qmin = -0.01;      % global min and max values for data
    %              qmax = 1.01;
    %
    %              q(q == qmax) = qmax*(1-1e-8);     % So floor works below.
    %              slope = (q - qmin)./(qmax-qmin);  % slopes must be in [0,1]
    %
    %              idx(m1) = floor(1 + slope(m1)*(n1 - 1));
    %              idx(m2) = n1 + floor(1 + slope(m2)*(n2 - 1));
    %
    %              set(p,'cdata',idx);
    %              fv = get(p,'FaceVertexCData');
    %
    %              if (sum(isnan(fv)) > 0)
    %                error('setcolors : Not all values in q have been mapped')
    %              end;
    %
    %              set(p,'FaceVertexCData',cm_new(fv,:));
    %              set(p,'facecolor','flat');
    %
    %
    %      See also PATCH, COLORMAP, UNDEROVER, UNDEROVER_COLORBAR,
    %      MULTICOLORMAP, MULTICOLORMAP_COLORBAR.
    %
    
    % cm = colormap;
    Ncm = 256;
    fine_pos = linspace(0, 1, Ncm)';
    idx = 0*q + nan;

    % m1 is dry land
    mt = isnan(q);
    mq = ~mt;

    % Load DEM data
    % filename = 'Missoula4Brian/topo/Scenario4a_maxice_except_Okanogan/topo_with_ice.tt3'; 
    filename = 'results_2_5_2a/resampled.asc'

    % plot topography
    [Z,tmin,tmax] = plot_topo(x,y,filename);

    % q colormap
    mymap = [ 127 255 0 
    0 128 0 
    255 255 0
    255 165 9
    255 0 0
    255 192 203
    255 0 255]/255;

    pos = linspace(0, 1, size(mymap, 1));
    
    % Create a q-grained colormap using interpolation
    
    q_mymap = interp1(pos, mymap, fine_pos);

    % q_mymap = createcolormap(Ncm,mymap(1,:),mymap(2,:),mymap(3,:),mymap(4,:),mymap(5,:), ...
    %     mymap(6,:),mymap(7,:));
   
    qmin = 0.00;
    qmax = 1295;
    
    q(q < qmin) = qmin;
    q(q > qmax) = qmax;
    
    slope = (q - qmin)./(qmax-qmin);
    
    idx(mq) = round(Ncm +1 + slope(mq)*(Ncm - 1));
    
    
    % Map nan values to brown (land)
    % cnew = [[153,76,0]/255; cm];
    terrain = [102 205 170
        255 255 0
        0 128 0
        255 165 0
        139 0 0
        165 42 42
        128 128 128
        255 250 250]/255;
    alpha = 0.0;
    terrain  = (1-alpha)*terrain + alpha;
    pos_t = linspace(0,1,size(terrain,1));
    t_map = interp1(pos_t,terrain,fine_pos);
  
    % t_map = createcolormap(Ncm,terrain(1,:),terrain(2,:),terrain(3,:),terrain(4,:),terrain(5,:), ...
    %     terrain(6,:),terrain(7,:),terrain(8,:)); % 256x3 array

    tmin = -115;
    tmax = 3731;

    Z(Z < tmin) = tmin;

    Z(Z > tmax) = tmax;

    slope = (Z - tmin)./(tmax-tmin);
    
    idx(mt) = round(1 + slope(mt)*(Ncm - 1));

    % done with mapping indices
    cm = [t_map;q_mymap];

    % cnew = cm;
    colormap(cm);
    
    set(p,'cdata',idx);
    fv = get(p,'FaceVertexCData');
    
    if (sum(isnan(fv)) > 0)
        fprintf('setcolors : Not all values in q have been mapped\n')
    end
    
    if (min(fv(:)) < 1)
        fprintf('Error in min index : \n');
    end
    
    
    if (max(fv(:)) > length(cm))
        fprintf('Error in max index : \n');
    end
    
    set(p,'FaceVertexCData',cm(fv,:));
    set(p,'facecolor','flat');
    
    set(p,'cdatamapping','direct');

    return;
    
end


function [Z,tmin,tmax] = plot_topo(xp,yp,filename)

    % Read the header information
    fid = fopen(filename, 'r');

    ncols = fscanf(fid,'%d',1); fscanf(fid,'%s',1);
    nrows = fscanf(fid,'%d',1); fscanf(fid,'%s',1);
    xllcorner = fscanf(fid,'%g',1);   fscanf(fid,'%s',1);
    yllcorner = fscanf(fid,'%g',1);   fscanf(fid,'%s',1);
    cellsize = fscanf(fid,'%g',1);    fscanf(fid,'%s',1);
    NODATA_value = fscanf(fid,'%g',1); fscanf(fid,'%s',1);
    T = fscanf(fid,'%g',nrows*ncols);
    fclose(fid);

    % Read the elevation data
    T = reshape(T,ncols,nrows);
    elevation = fliplr(T);
    
    % Create x and y coordinates for each elevation point
    x = xllcorner + (0:ncols-1) * cellsize;
    y = yllcorner + (0:nrows-1) * cellsize;

    % Create a grid for the elevation data
    [X, Y] = meshgrid(x, y);
    
    tmin = min(elevation(:));
    tmax = max(elevation(:));
    Z = interp2(X,Y,elevation',xp,yp,"linear");

end