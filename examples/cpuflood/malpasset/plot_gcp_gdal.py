# Title: Using GCPs to georeference a raster using GDAL
# Author: Brian KYANJO
# Date: 02/08/2023 10:00:00 am

# === Import libraries ===
import pandas as pd
import numpy as np
from osgeo import gdal
import os

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

def rewrite_kml(kml_file):
    """
    This function rewrites the *.kml file to change the latlon box, range, and average lat and lon
    """
    #computational domain (latlon box) for the image
    north = '1848572.75'
    south = '1832407.25'
    east = '959554.0'
    west = '953236.0'

    # replace with the new values(different gcp points have different latlon box)
    north_r = '43.548787426'
    south_r = '43.400118651'
    east_r = '6.781588849'
    west_r = '6.689219689'

    # compute the average lat and lon
    av_lat = (float(north) + float(south))/2
    av_lon = (float(east) + float(west))/2
    
    average_lat = (float(north_r) + float(south_r))/2
    average_lon = (float(east_r) + float(west_r))/2

    #compute the range of lat and lon in meters
    range = range_kml(south,north,west,east) 
    range_r  = range_kml(south_r,north_r,west_r,east_r)

    # replace the range
    with open(kml_file, 'r') as f:
        filedata = f.read()
        file = filedata.replace(str(range), str(range_r))
    
    # write the new kml file
    with open(kml_file, 'w') as f:
        f.write(file)

    # replace the averages
    with open(kml_file, 'r') as f:
        filedata = f.read()
        file = filedata.replace(str(av_lat), str(average_lat))

    with open(kml_file, 'w') as f:
        f.write(file)   

    with open(kml_file, 'r') as f:
        filedata = f.read()
        file = filedata.replace(str(av_lon), str(average_lon))

    with open(kml_file, 'w') as f:
        f.write(file)

    # replace the latlon box
    with open(kml_file, 'r') as f:
        filedata = f.read()
        file = filedata.replace(north, north_r)

    # write the new kml file
    with open(kml_file, 'w') as f:
        f.write(file)

    with open(kml_file, 'r') as f:
        filedata = f.read()
        file = file.replace(south, south_r)

    # write the new kml file
    with open(kml_file, 'w') as f:
        f.write(file)

    with open(kml_file, 'r') as f:
        filedata = f.read()
        file = file.replace(east, east_r)

    # write the new kml file
    with open(kml_file, 'w') as f:
        f.write(file)

    with open(kml_file, 'r') as f:
        filedata = f.read()
        file = file.replace(west, west_r)

    # write the new kml file
    with open(kml_file, 'w') as f:
        f.write(file)

# check if the _plots folder exists
if os.path.exists("_plots"):
    print("_plots folder exists")
else:
    print("_plots folder does not exist")
    os.system("python make_plots_kml.py")

# go into the _plots folder
os.chdir("_plots")

# check if the kmz folder exists
if os.path.exists("kmz"):
    print("kmz folder exists")
else:
    print("kmz folder does not exist")
    os.system("mkdir kmz")

os.chdir("kmz") # go into the kmz folder
os.system("cp -f ../*.kmz .") # copy the kmz file to the kmz folder
os.system("unzip -o *.kmz") # unzip the kmz file, overwrite if unzipped already (to avoid georeferencing twice)

os.chdir("fig1") # go into the fig1 folder

# locate the gcp points file
gcp_points = "../../../../malpasset_gcp.points" # locate the gcp points file from the frame folders

print('Geo-referencing images...')

# access images in the frame* folders
for folder in os.listdir():
    if folder.startswith("frame"):
        os.chdir(folder) # go into the frame* folder
        rewrite_kml("doc.kml") # edit the *.kml file
        for image in os.listdir():
            if image.endswith("png"):
                image_in = image # get the image name
                image_out = gcp_gdal(gcp_points, image_in) # call the gcp_gdal function
                os.system("rm " + image_in[:-4] + "_gcp.tif") # remove the _gcp image (size is large)
        os.chdir("..") # go back to the fig1 folder

# go back to the kmz folder
os.chdir("..")
rewrite_kml("doc.kml") # edit the *.kml file

# finally open the .kml file to visualize the georeferenced image in google earth
print("Opening the .kml file in Google Earth")
os.system("open doc.kml")