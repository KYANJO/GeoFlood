import numpy as np

def read_topo_data(topofile):

    f = open(topofile,'r')

    l = f.readline()
    ncols = np.fromstring(l.split()[0].strip(),sep=' ')                    

    l = f.readline()
    nrows = np.fromstring(l.split()[0].strip(),sep=' ')

    l = f.readline()
    xllcorner = np.fromstring(l.split()[0].strip(),sep=' ')

    l = f.readline()
    yllcorner = np.fromstring(l.split()[0].strip(),sep=' ')

    l = f.readline()
    cellsize = np.fromstring(l.split()[0].strip(),sep=' ')

    return ncols[0],nrows[0],xllcorner[0],yllcorner[0],cellsize[0]

def compute_distances(lon,lat):

    # lon = [ax bx bx ax ax];
    # lat = [ay ay by by ay];

    lon = np.array([lon[0],lon[1],lon[1],lon[0],lon[0]])
    lat = np.array([lat[0],lat[0],lat[1],lat[1],lat[0]])

    labels = ['top', 'right','bottom','left']

    d = np.zeros(4)
    for i in range(0,4):
        lat1 = np.deg2rad(lat[i])
        lat2 = np.deg2rad(lat[i+1])
        lon1 = np.deg2rad(lon[i])
        lon2 = np.deg2rad(lon[i+1])
    
        dlon = lon2-lon1
        dlat = lat1-lat2
    
        R = 6367.5e3
        a = (np.sin(dlat/2.0))**2 + np.cos(lat1) * np.cos(lat2) * (np.sin(dlon/2.0))**2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
        d[i] = R * c

        print("%20s %12.4f" % (labels[i],d[i]))

    print("")
    avg_width = (d[0] + d[2])/2.0
    avg_height = (d[1] + d[3])/2.0
    print("%20s %12.4f" % ("Average width",avg_width))
    print("%20s %12.4f" % ("Average height",avg_height))

    print("")
    print("%20s %12.6f" % ("Ratio edge/avg width",d[1]/avg_width))
    print("%20s %12.6f" % ("Ratio top/bottom",d[0]/d[2]))
    print("%20s %12.6f" % ("Ratio left/right",d[1]/d[3]))

    return d


def deg2m2(lon,lat):

    area_deg = (lat[1]-lat[0])*(lon[1]-lon[0])
    d = compute_distances(lon,lat)

    avg_width = (d[0] + d[2])/2.0
    avg_height = (d[1] + d[3])/2.0

    area_m2 = avg_width*avg_height

    return area_m2/area_deg


def region_coords(xll,xur, num_cells,lower,upper):
    """ Get location of region in coarse grid coordinates
    xll, xur     : lower and upper coordinates of approximate region
    num_clels    : Number of grid cells at coarsest level
    lower, upper :
    """

    # Get degrees per finest level :
    dx_coarse = (upper[0]-lower[0])/num_cells[0]
    dy_coarse = (upper[1]-lower[1])/num_cells[1]

    # Zoom region
    # xll = [-111.623926, 43.913661]
    # xur = [-111.620150, 43.916382]
    mx_xlow = np.floor((xll[0] - lower[0])/dx_coarse).astype(int)
    my_ylow = np.floor((xll[1] - lower[1])/dy_coarse).astype(int)
    xlow = lower[0] + mx_xlow*dx_coarse
    ylow = lower[1] + my_ylow*dy_coarse

    mx_zoom = np.ceil((xur[0] - xlow)/dx_coarse).astype(int)
    my_zoom = np.ceil((xur[1] - ylow)/dy_coarse).astype(int)

    region_lower = [xlow,ylow]
    region_upper = [xlow + mx_zoom*dx_coarse, ylow + my_zoom*dy_coarse]

    figsize = np.array([mx_zoom, my_zoom])   # [1,1]

    return region_lower, region_upper, figsize

def read_locations_data(malpasset_loc):
    with open (malpasset_loc, "r") as myfile:
        data = myfile.read().splitlines()
    data2 = []
    for i in range(1,len(data)):
        data[i] = data[i].split()
        data2.append([float(j) for j in data[i]])

    data2 = np.array(data2)
    x = []
    y = []
    for data in data2:
        x.append(data[1])
        y.append(data[2])
    police = [range(17),x[:17], y[:17]]
    transformers = [range(17,20),x[17:20], y[17:20]]
    gauges = [range(6,15),x[20:29], y[20:29]]
    all_guages = [range(29),x, y]
    
    return police, transformers, gauges, all_guages 

def ozone(x2,y2,z2):
    """
    @description: Convert ECEF coordinates to geodetic coordinates using Ozone's method
    @source: obtained from https://www.researchgate.net/publication/254250446_Transforming_cartesian_coordinates_X_Y_Z_to_geographical_coordinates_h 
    """

    a = 6378137 # WGS-84 ellipsiod semi-major axis
    f = 1/298.257223563 # WGS-84 ellipsiod flattening parameter
    b = a*(1-f) # WGS-84 ellipsiod semi-minor axis
    p = np.sqrt((x2**2) + (y2**2)) # distance from z-axis
    e = np.sqrt(f*(2 - f)) # eccentricity
    #==== Ozone's method =====
    M = (a*p - ((a**2) - (b**2)))/(2*b*z2)
    N = (a*p + ((a**2) - (b**2)))/(2*b*z2)
    V = 4*N*M + 1
    W = 2*(N**2 - M**2)
    I = (np.sqrt((V/3)**3 + (W/2)**2) + (W/2))**(1/3) - (np.sqrt((V/3)**3 + (W/2)**2) - (W/2))**(1/3)
    J = np.sqrt(2*I + 4*M**2)
    K = (2*(N - M*I))/J
    G = (2*M + J)**2 - 4*(I - K)
    u = (2*M + J + np.sqrt(G))/2
    t_phi = 2*a*u/(b*(u**2 - 1))
    lat = np.degrees(np.arctan(t_phi))
    lon = np.degrees(np.arctan(y2/x2))
    v = a/(np.sqrt(1 - (e**2)*np.degrees(np.sin(lat)**2))) # radius of curvature in the prime vertical
    h = (p/np.degrees(np.cos(lat))) - v
    return lat, lon, h

def xyz_to_latlon(file):   
    with open (file, "r") as myfile:
        data = myfile.read().splitlines()
    data2 = []
    for i in range(2, len(data)):
        data[i] = data[i].split()
        data2.append([float(j) for j in data[i]])

    data2 = np.array(data2)
    x = []; y = []; z = []
    for data in data2:
        x.append(data[1])
        y.append(data[2])
        z.append(data[3])

    R = 6371e3 # metres radius of earth

    #====== X,Y,Z are obtained from https://tool-online.com/en/coordinate-converter.php took the corresponding XYZ for latlon coordinates of the reservior ======
    # used [lat,lon,h] = [43.512128, 6.75660, 55m] of the reservior as the starting point and converted it to XYZ using the above website to obtain the following XYZ of the reservior

    [X,Y,Z] = [4600827.391, 545086.101,  4368976.942] #XYZ of the reservior

    x2 = [X]; y2 = [Y]; z2 = [Z]
    for i in range(1,len(x)):
        x1 = x[i] - x[i-1] #distance between two x points
        y1 = y[i] - y[i-1] #distance between two y points
        z1 = z[i] - z[i-1] #distance between two z points
        X = X + x1
        Y = Y + y1
        Z = Z + z1
        x2.append(X)
        y2.append(Y)
        z2.append(Z)

    x2 = np.floor(x2).astype(int)
    y2 = np.floor(y2).astype(int)
    z2 = np.floor(z2).astype(int)

    #===== converting xyz to latlon ======
    lat, lon, h = ozone(x2,y2,z2)

    #===== write to a file ======
    with open('latlon.csv', 'w') as f:
        for i in range(len(lat)-1):
            f.write(str(lat[i])+ " " + str(lon[i]) + " " + str(h[i]))
            f.write("\n")

    return lat, lon, h

