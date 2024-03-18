# @author :: Brian Kyanjo
# @data   :: January 21, 2024

import numpy as np
import os

def read_topo_data(topofile):

    f = open(topofile,'r')
    l = f.readline()
    first_entry = l.split()[0]
    
    if first_entry.isdigit():
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
    else:
        ncols = np.fromstring(l.split()[1].strip(),sep=' ')

        l = f.readline()
        nrows = np.fromstring(l.split()[1].strip(),sep=' ')

        l = f.readline()
        xllcorner = np.fromstring(l.split()[1].strip(),sep=' ')

        l = f.readline()
        yllcorner = np.fromstring(l.split()[1].strip(),sep=' ')

        l = f.readline()
        cellsize = np.fromstring(l.split()[1].strip(),sep=' ')

        return ncols[0],nrows[0],xllcorner[0],yllcorner[0],cellsize[0]


# =================== swap headers ======================================
def swap_header(fp):

    f    = open(fp, 'r')
    l    = f.readline()
    firstline = l.split()[0]

    if firstline.isdigit(): # check if firstline is a number
        
        fnew = open(fp + '_new', 'w')

        # swap the header use a new name for the entry
        fnew.write("ncols %s\n" % (l.split()[0]))

        l = f.readline()
        fnew.write("nrows %s\n" % (l.split()[0]))

        l = f.readline()
        fnew.write("xllcorner %s\n" % (l.split()[0]))

        l = f.readline()
        fnew.write("yllcorner %s\n" % (l.split()[0]))

        l = f.readline()
        fnew.write("cellsize %s\n" % (l.split()[0]))

        l = f.readline()
        fnew.write("NODATA_value %s\n" % (l.split()[0]))

        # write the rest of the file
        for line in f:
            fnew.write(line)

        fnew.close()
        f.close()
        return fnew

    else:
        fnew = open(fp + '_new', 'w')
        #  rename headers regardless
        fnew.write("ncols %s\n" % (l.split()[1]))
        l = f.readline()
        fnew.write("nrows %s\n" % (l.split()[1]))
        l = f.readline()
        fnew.write("xllcorner %s\n" % (l.split()[1]))
        l = f.readline()
        fnew.write("yllcorner %s\n" % (l.split()[1]))
        l = f.readline()
        fnew.write("cellsize %s\n" % (l.split()[1]))
        l = f.readline()
        fnew.write("NODATA_value %s\n" % (l.split()[1]))

        for line in f:
            fnew.write(line)
        fnew.close()
        f.close()
        return fnew

# =================== merge dems ======================================
def merge_dems(dem_paths, output_filename):
    """
    Merge multiple DEMs into a single file.

    Parameters
    ----------
    dem_paths : list
        A list of file paths to the DEMs to be merged.
    output_filename : str
        The file path to the output file.

    Returns
    -------
    None
    """

    # import the necessary modules
    import rasterio
    from rasterio.merge import merge

    # merge the DEMs
    src_files_to_mosaic = []
    for fp in dem_paths:
        swap_header(fp) # swap the header of the DEM if necessary
        src = rasterio.open(fp + '_new')   # open the new DEM
        src_files_to_mosaic.append(src)   # append the DEM to the list of DEMs to be merged
        mosaic, out = merge(src_files_to_mosaic)  # merge the DEMs

    # copy the file meta data to the mosaic
    out_meta = src.meta.copy()
    out_meta.update(
        {'driver': 'GTiff',
            'height': mosaic.shape[1],
            'width': mosaic.shape[2],
            'transform': out,
            }
    )

    # write the mosaic to a new GeoTIFF file
    tif_file = output_filename[:-4] + '.tif'
    with rasterio.open(
            tif_file,
            "w",
            **out_meta
    ) as dest:
        dest.write(mosaic)

    # ==== save the mosaic as a new ASCII DEM file ========================
    # Construct the gdal_translate command
    command = [
        'gdal_translate',
        '-of', 'AAIGrid',  # Output format
        tif_file,  # Input TIFF file
        output_filename  # Output ASCII DEM file
    ]

    # Execute the command using subprocess
    import subprocess
    subprocess.run(command, check=True)

    # remove the temporary files with errors
    for fp in dem_paths:
        if os.path.isfile(fp + '_new'): os.remove(fp + '_new')
        if os.path.isfile(fp[:-4] + '_merged.tif'): os.remove(fp[:-4] + '_merged.tif')
        if os.path.isfile(fp + '.xml'): os.remove(fp + '.xml')
        if os.path.isfile(fp[:-4] + '_merged.tt2.aux.xml'): os.remove(fp[:-4] + '_merged.tt2.aux.xml')
            


def convert_file_type(input_file,output_file,input_type,output_type):
    """
    Reads a file and converts from one type to another
        file_type = 1 : columns of data: x y z
        file_type = 2 and 3 : header describes grid format followed by z data 
            file_type = 2 : one value per line
            file_type = 3 : mx values per line
    returns a converted file
    """

    # check if output file exists
    if os.path.isfile(output_file):
        print(f'Warning: Converted file: {output_file} already exists, skipping file conversion, delete it to overwrite.')
        return
    else:
        print (f'Converting {input_file} from type {input_type} to type {output_type} and saving as {output_file}.')

    # check input type
    if input_type == output_type:
        return input_file  
    elif input_type == 2 or input_type == 3:
        # read in header information
        mx,my,xll,yll,cellsize = read_topo_data(input_file)
        mx = int(mx)
        my = int(my)
        dx = cellsize
        dy = dx             # grid spacing  

        # read in the remaining data
        if input_type == 2:
            data = np.loadtxt(input_file,skiprows=6)
            # save data as 1D array
            data = data.flatten()
        
        elif input_type == 3:
            # convert 3 to 2
            f = open(input_file,'r')
        
            # skip six header lines
            for i in range(6):
                f.readline()

            data = [] # save z values as a list
            for j in range(my):
                l = f.readline()
                data2str = l.split()
                for k in range(mx):
                    data.append(float(data2str[k]))

    elif input_type == 1:
        print('Error: Converting from file type 1 to file type 2 or 3 is not yet supported.')

    
    # convert from 2 or 3 to 1
    if input_type == 2 or input_type == 3 and output_type == 1:
        # write out the data
        f = open(output_file,'w')
        k = 0
        for j in range(my):
            for i in range(mx):
                # f.write('%f %f %f\n' % (xll+i*dx,yll+j*dy,data[j*mx+i]))
                f.write('%f %f %f\n' % (xll+i*dx,yll+j*dy,data[k]))
                k = k + 1
        f.close()
    elif input_type == 3 and output_type == 2:
        # write out the data
        f = open(output_file,'w')
        f.write('%d\tncols\n' % mx)
        f.write('%d\tnrows\n' % my)
        f.write('%f\txllcorner\n' % xll)
        f.write('%f\tyllcorner\n' % yll)
        f.write('%f\tcellsize\n' % dx)
        f.write('%f\tnodata_value\n' % -9999.0)
        for idx in range(my*mx):
            # for i in range(mx):
            f.write('%f\n' % data[idx])
        f.close()
               
    return output_file

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


