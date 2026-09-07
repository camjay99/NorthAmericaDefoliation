import argparse
import json

import ee

import utils.geometries as geometries
import utils.preprocessing as preprocessing
import utils.submission as submission


##############################################################
# Parse arguments
##############################################################

parser = argparse.ArgumentParser(
    description='Options for calculating growing season')

# The script will ONLY submit the run when -s or --submit is included.
parser.add_argument('--submit', '-s', action='store_true')

# Cloud storage bucket to save results in.
parser.add_argument('--bucket', '-b', action='store', default=None)

# The project to submit the code in. 
# You may be prompted to to authenticate.
parser.add_argument('--project', '-p', action='store', 
                    default=None, required=True)

# The crs to use for the output
parser.add_argument('--crs', '-c', action='store', default='epsg:5070')

# The width/length of grid cells to use for computation (in lat/lon degrees)
parser.add_argument('--width', '-w', action='store', type=float, default=0.75)
parser.add_argument('--length', '-l', action='store', type=float, default=0.75)


# Parse arguments provided to script
args = parser.parse_args()


##############################################################
# Initialize Google Earth Engine API
##############################################################

try:
    ee.Initialize(project=args.project)
except:
    # need to authenticate with your credential at the first time
    ee.Authenticate()
    ee.Initialize(project=args.project)


##################################################################
# Specify base names and load previous results
##################################################################

assert (args.bucket is not None), "Must specify bucket if exporting to cloud storage."
file_name_prefix = f'topography_North_America/topography'
image_manifests = {}
description_base = f'North_America_Topography'
geometry = geometries.get_range()


##################################################################
# Split study regions into grid cells of specified size.
##################################################################

#Specify grid size in projection, x and y units (based on projection).
projection = 'EPSG:4326' # WGS84 lat lon

# Make grid and visualize.
proj = ee.Projection(projection).scale(args.width, args.length)
grid = geometry.coveringGrid(proj)

gridSize = grid.size().getInfo()
gridList = grid.toList(gridSize)

for i in range(gridSize):
    gridCell = ee.Feature(gridList.get(i)).geometry()

    ##################################################################
    # Prepare Data
    ##################################################################

    # Topography
    elevation = ee.Image('USGS/SRTMGL1_003')
    slope = ee.Terrain.slope(elevation).rename('slope')
    aspect = ee.Terrain.aspect(elevation).rename('aspect')

    # Landcover
    landcover = ee.Image("USGS/NLCD_RELEASES/2020_REL/NALCMS")
    forest = landcover.gte(41).And(landcover.lte(43)).rename('forest')
    ## Distance to forest edge
    kernel = ee.Kernel.euclidean(radius=4000, units='meters', normalize=False)
    distance_to_forest = (forest.Not()
                          .distance(kernel)
                          .rename('distance_to_edge'))

    output = ee.Image([elevation, slope, aspect, landcover, distance_to_forest])

    ## TODO: Compress output data to more reasonable size.


    ##################################################
    # Export results
    ##################################################

    if args.submit:
        submission.submit_job(
            image=output,
            assetID='',
            file_name_prefix=file_name_prefix,
            description_base=description_base,
            scale=preprocessing.resolutions[args.data],
            crs=args.crs,
            region=gridCell,
            cloudstorage=args.cloudstorage,
            bucket=args.bucket,
            i=i
        )