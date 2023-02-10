# Title: Using GCPs to georeference a raster using GDAL
# Author: Brian KYANJO
# Date: 02/08/2023 10:00:00 am

# === Import libraries ===
import os
import numpy as np
import pandas as pd
from time import sleep
from osgeo import gdal
from random import random
from multiprocessing import Pool

# === Define functions ===
def gcp_gdal(gcp_points, image_in):
   
    """
    This function reads in a gcp file and an image and creates a georeferenced image.
    check this website  https://kokoalberti.com/articles/georeferencing-and-digitizing-old-maps-with-gdal/ for details about using gdal to georeference images
    """
   
    # === read in gcp points file ===
    data = pd.read_csv(gcp_points)  # read in gcp points file
    data = np.array(data)           # convert the data to a numpy array
    data = data[1:,:8]              # extract wanted columns and rows
    data = pd.DataFrame(data)       # convert the data to a dataframe
    data.columns = ['mapX','mapY','pixelX', 'pixelY', 'enable', 'dX','dY', 'residual'] # rename the columns
    data = data.drop(columns=['enable', 'dX','dY', 'residual']) # drop unwanted columns
    data['pixelY'] = data['pixelY'].astype(float)*(-1)          # flip the y axis otherwise the image will be upside down

    # === create a list of gcp points ===
    gcp_list = [] # create an empty list

    for index, rows in data.iterrows():
        gcp_list.append( "-gcp" +" "+ rows.pixelX +" "+ str(rows.pixelY) +" "+ (rows.mapX) +" "+ (rows.mapY)) # append the gcp points to the list

    # === create a string of gcp points ===
    gcp_ = " "                       # create an empty string
    for i in range(len(gcp_list)):
        gcp_ += (gcp_list[i]) + " "  # append the gcp points to the string

    # === create a georeferenced image using GDAL ===
    image_gcp = image_in[:-4] + "_gcp.tif" # create a new image name
    os.system("gdal_translate" + gcp_ + "-of GTiff "+" "+ image_in +" "+ image_gcp) # create a new GeoTiff image using gdal_translate

    os.system("gdalinfo" + " " + image_gcp) # check for the GCPs now in our new file with gdalinfo

    image_projected = image_in[:-4] + ".tif" # create a new image name
    os.system("gdalwarp -tps -r near -t_srs EPSG:4326 -overwrite -co COMPRESS=JPEG -co JPEG_QUALITY=50 -co TILED=YES" +" "+ image_gcp +" "+ image_projected) # reprojecting the image to EPSG:4326, compressing, and Tiling it using gdalwarp

    os.system("gdaladdo --config COMPRESS_OVERVIEW JPEG --config PHOTOMETRIC_OVERVIEW YCBCR  --config JPEG_QUALITY_OVERVIEW 50 -r average " +" " + image_projected +" "+ "2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536") # create overviews for the image using gdaladdo (low resolution images)

    image_out = image_in[:-4] + ".png" # create a new image name
    os.system("mv " +" "+ image_projected +" "+ image_out) # rename the image

    return image_out # return the georeferenced image

def range_kml(y1,y2,x1,x2):
    '''(Range is the Distance in meters from the point specified by <longitude>, <latitude>, and <altitude> to the LookAt position.)
    '''
    R = 6371 # radius of earth
    domain_width = R*np.cos(abs(float(y1)+float(y2))*np.pi/360.0)*(float(x2)-float(x1))*np.pi/180.0
    dist_factor = 2  # factor by which height should exceed width
    initial_height = min([1000*dist_factor*domain_width,9656064.0])   # <= 6000 miles
    return initial_height

def read_locations_data(malpasset_loc):
    """
    This function reads in the locations of the police stations, transformers, and gauges from a file
    """
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
    
    return police, transformers, gauges

def rewrite_kml(kml_file,coordinates,malpasset_loc,gauge_lat_long):
    """
    This function rewrites the *.kml file to change the latlon box, range, and average lat and lon
    """
    #computational domain (latlon box) for the image
    [north,south,east,west] = coordinates[0]

    # in exponential and interger notation
    north_e = '{:.5e}'.format(float(north))
    south_e = '{:.5e}'.format(float(south))
    east_e = str(int(float(east)))
    west_e = str(int(float(west)))

    # replace with the new values(different gcp points have different latlon box)
    [north_r, south_r, east_r, west_r] = coordinates[1]

    # since -180 < east_r, west_r > 180 then
    east_i =  "{:.4f}".format(float(str(float(east) - 360)))
    west_i = "{:.4f}".format(float(str(float(west) - 360)))

    # compute the average lat and lon
    av_lat = (float(north) + float(south))/2
    av_lon = (float(east) + float(west))/2
    
    average_lat = (float(north_r) + float(south_r))/2
    average_lon = (float(east_r) + float(west_r))/2

    #compute the range of lat and lon in meters
    range_ = range_kml(south,north,west,east) 
    range_r  = range_kml(south_r,north_r,west_r,east_r)

    # check if in the kml directory 
    if os.path.exists("levels"): # choose levels since its the only directory in the kml directory

        # approximate the gauge locations
        [gauge_lat,gauge_lon] = gauge_lat_long

        police, transformers, gauge = read_locations_data(malpasset_loc)

        gauge_loc = []
        gauge_round = []
        for i in range(len(gauge[0])):
            gauge_loc.append([gauge[1][i],gauge[2][i]])  
            y2 = "{:.5e}".format(gauge_loc[i][1])
            if y2[6] == '0':
                y2 = "{:.4e}".format(gauge_loc[i][1])
            gauge_round.append([int(gauge_loc[i][0].round()),y2])
        
        #for x1,x2 > 180 and < -180
        gauge_360 = np.array(gauge_loc) - 360
        
        # gauge list
        gauge_list = [gauge_loc,gauge_round,gauge_360]
        # replace the gauge locations
        for i in range(len(gauge_loc)):
            for j in range(len(gauge_list)):
                with open(kml_file, 'r') as f:
                    filedata = f.read()
                    file = filedata.replace(str(gauge_list[j][i][0]), str(gauge_lon[i]))
                with open(kml_file, 'w') as f:
                    f.write(file)

                with open(kml_file, 'r') as f:
                    filedata = f.read()
                    file = filedata.replace(str(gauge_list[j][i][1]), str(gauge_lat[i]))
                with open(kml_file, 'w') as f:
                    f.write(file)

    # replace the range and lat-long averages
    av_lat_s = [str(range_),str(av_lat),str(av_lon)]
    av_lat_r = [str(range_r), str(average_lat),str(average_lon)]
    for i in range(len(av_lat_s)):
        with open(kml_file, 'r') as f:
            filedata = f.read()
            file = filedata.replace(av_lat_s[i], av_lat_r[i])
        with open(kml_file, 'w') as f:
            f.write(file)

    # replace the latlon box
    domain_box = [north,south,east,west]
    domain_box_i = [north,south,east_i,west_i] #for simpilicity
    domain_boxe = [north_e,south_e,east_e,west_e]
    domain_search = [domain_box,domain_box_i,domain_boxe]
    domain_box_replace = [north_r,south_r,east_r,west_r]
    for i in range(len(domain_search)):
        for j in range(len(domain_box_replace)):
            with open(kml_file, 'r') as f:
                filedata = f.read()
                file = filedata.replace(domain_search[i][j], domain_box_replace[j])
            with open(kml_file, 'w') as f:
                f.write(file)

def overlay_image_google_earth(func_arg):
    '''
    This function overlays the a georeferenced simulated image onto google earth.
    '''
    # get the arguments
    gcp_points,coordinates,malpasset_loc,gauge_lat_long = func_arg

    # check if the _plots folder exists
    if not os.path.exists("_plots"):
        print("_plots folder does not exist, creating _plots folder...")
        os.system("python make_plots_kml.py")

    # go into the _plots folder
    os.chdir("_plots")

    # check if the kmz folder exists
    if not os.path.exists("kmz"):
        print("kmz folder does not exist, creating kmz folder...")
        os.system("mkdir kmz")

    os.chdir("kmz") # go into the kmz folder
    os.system("cp -f ../*.kmz .") # copy the kmz file to the kmz folder
    os.system("unzip -o *.kmz") # unzip the kmz file, overwrite if unzipped already (to avoid georeferencing twice)

    os.chdir("fig1") # go into the fig1 folder

    print('Geo-referencing images...')

    # access images in the frame* folders
    for folder in os.listdir():
        if folder.startswith("frame"):
            os.chdir(folder) # go into the frame* folder
            rewrite_kml("doc.kml",coordinates,malpasset_loc,gauge_lat_long) # edit the *.kml file
            for image in os.listdir():
                if image.endswith("png"):
                    image_in = image # get the image name
                    image_out = gcp_gdal(gcp_points, image_in) # call the gcp_gdal function
                    os.system("rm " + image_in[:-4] + "_gcp.tif") # remove the _gcp image (size is large)
            os.chdir("..") # go back to the fig1 folder

    os.chdir("../kml")  # go back to the kml folder
    rewrite_kml("regions.kml",coordinates,malpasset_loc,gauge_lat_long) # rewrite the regions.kml file
    rewrite_kml("gauges.kml",coordinates,malpasset_loc,gauge_lat_long) # rewrite the gauges.kml file
    os.chdir("..") # go back to the kmz folder
    rewrite_kml("doc.kml",coordinates,malpasset_loc,gauge_lat_long) # edit the *.kml file

    # make a .kmz file
    print("Making a .kmz file")
    os.system("rm -f *.kmz") # remove the old kmz file

    # check if the _geoReferenced folder exists
    if not os.path.exists("../_geoReferenced"):
        os.mkdir("../_geoReferenced")
   
    os.system("zip -r -o -9 ../_geoReferenced/MalpassetDam.kmz *") # make a new .kmz file and save it in the _geoReferenced folder
    os.chdir("../_geoReferenced") # go back to the _geoReferenced folder

    # finally open the .kmz file to visualize the georeferenced image in google earth
    print("Opening the .kmz file in Google Earth")
    os.system("open MalpassetDam.kmz")

# === end of function definitions ===

# === define the latlon box ===
north = '1848572.75'
south = '1832407.25'
east = '959554.0'
west = '953236.0'

# === define the latlon box for the reference image ===
north_r = '43.549942109'
south_r = '43.400191065'
east_r = '6.781941194'
west_r = '6.690324187'

# === guages locations latlong (approximate) ===
gauge_lat = [43.508383,43.503714,43.493864,43.488525,43.476069,43.466011,43.450400,43.432585,43.428329]
gauge_lon = [6.757204,6.758799,6.751142,6.740853,6.743523,6.736402,6.735097,6.721494,6.716970]

# === function arguments ===
gauge_lat_long = [gauge_lat, gauge_lon]   
coordinates = [[north,south,east,west],[north_r,south_r,east_r,west_r]]    

# === locate data files ===
malpasset_loc = "../../../malpasset_locs.txt"   # Police, transformer and guage data
gcp_points = "../../../../malpasset_gcp.points" # locate the gcp points file from the frame folders

# make one fuction argument
func_arg = [gcp_points,coordinates,malpasset_loc,gauge_lat_long]

def task(func_arg):
    # unpack the function arguments
    gcp_points,coordinates,malpasset_loc,gauge_lat_long = func_arg
    # block for a moment
    sleep(random())
    # call the function
    return overlay_image_google_earth(func_arg)

# === parallelization using multiprocessing ===
if __name__ == '__main__':
    with Pool(processes = 4) as p:
        # prepare the arguments
        args = [func_arg]
        # run the function
        p.map(task, args)


