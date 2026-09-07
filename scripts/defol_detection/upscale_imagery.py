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
    description='Options for calculating defoliation')

# The script will ONLY submit the run when -s or --submit is included.
parser.add_argument('--submit', '-s', action='store_true')

# Whether to export results to a cloud storage bucket. If true,
# `bucket` must also be set.
parser.add_argument('--cloudstorage', '-C', action='store_true')

# Cloud storage bucket to save results in.
parser.add_argument('--bucket', '-b', action='store', default=None)

# The project to submit the code in. 
# You may be prompted to to authenticate.
parser.add_argument('--project', '-p', action='store', 
                    default=None, required=True)

# The crs to use for the output
parser.add_argument('--crs', '-c', action='store', default='epsg:5070')

# The scale to use for the output
parser.add_argument('--scale', '-S', action='store', type=float, default=3000)

# Year to upscale imagery for.
parser.add_argument('--year', '-y', action='store', type=int, default=2021)

# The width/length of grid cells to use for computation (in lat/lon degrees)
parser.add_argument('--width', '-w', action='store', type=float, default=2)
parser.add_argument('--length', '-l', action='store', type=float, default=2)

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

if args.cloudstorage:
    assert (args.bucket is not None), "Must specify bucket if exporting to cloud storage."
    file_name_prefix = f'upscaled_imagery/upscaled_'
    image_manifests = {}
description_base = f'Upscaled'
assetID=f'projects/{args.project}/assets/upscaled_imagery/upscaled'


##################################################################
# Prepare imagery
##################################################################

defol_col = (ee.ImageCollection(f'projects/{args.project}/assets/defoliation_score_North_America')
             .filter(ee.Filter.eq('year', args.year))
             .filter(ee.Filter.eq('method', 'Theil-Sen SWIR1')));
nlcd_landcover = ee.ImageCollection('USGS/NLCD_RELEASES/2019_REL/NLCD') \
    .filter(ee.Filter.eq('system:index', '2019')).first().select('landcover')
esri_lulc_ts= (ee.ImageCollection("projects/sat-io/open-datasets/landcover/ESRI_Global-LULC_10m_TS")
               .filterDate('2019-01-01', '2020-01-01')
               .mosaic());

# Combine NLCD and ESRI LULC to ensure best maps where possible.
nlcd_valid = nlcd_landcover.gte(0).unmask(0, False) # base map for where NLCD is valid.
for_mask_us = nlcd_valid.And(nlcd_landcover.gte(41).And(nlcd_landcover.lte(43)).unmask(0, False))
for_mask_ca = nlcd_valid.Not().And(esri_lulc_ts.eq(2))
forest_mask = for_mask_us.Or(for_mask_ca)

forest_change = ee.Image("UMD/hansen/global_forest_change_2025_v1_13")
forest_change_mask = forest_change.select('lossyear').lte(args.year).unmask().Not()


##################################################################
# Prepare mesh for processing many image tiles
##################################################################

#Specify grid size in projection, x and y units (based on projection).
projection = 'EPSG:4326'; # WGS84 lat lon

# Make grid and visualize.
proj = ee.Projection(projection).scale(args.width, args.length)
grid = geometries.get_range().coveringGrid(proj)

gridSize = grid.size().getInfo()
gridList = grid.toList(gridSize)

for i in range(gridSize):
    gridCell = ee.Feature(gridList.get(i)).geometry()

    def simple_mask(image):
        defol_mask = image.neq(0)
        return (image.updateMask(forest_mask)
                     .updateMask(forest_change_mask)
                     .updateMask(defol_mask))
    
    defol_tile = (defol_col.filterBounds(gridCell)
                           .map(simple_mask)
                           .max()
                           .reproject('EPSG:5070', 
                                      [30, 0, 1130610, 0, -30, 2032620]))
    
    defol_tile = defol_tile.reduceResolution(
        ee.Reducer.mean(), 
        False, 
        15000)

    ##################################################
    # Export results
    ##################################################
    if args.submit:
        submission.submit_job(
            image=defol_tile,
            assetID=assetID,
            file_name_prefix=file_name_prefix,
            description_base=description_base,
            year=args.year,
            scale=args.scale,
            crs=args.crs,
            region=gridCell,
            cloudstorage=args.cloudstorage,
            bucket=args.bucket,
            i=i
        )
