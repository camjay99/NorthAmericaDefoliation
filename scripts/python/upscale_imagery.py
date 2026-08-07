import argparse
import json

import ee

import geometries
import preprocessing


##############################################################
# Parse arguments
##############################################################

parser = argparse.ArgumentParser(
    description='Options for calculating defoliation')

# The script will ONLY submit the run when -s or --submit is included.
parser.add_argument('--submit', '-s', action='store_true')

# The project to submit the code in. 
# You may be prompted to to authenticate.
parser.add_argument('--project', '-p', action='store', 
                    default=None, required=True)

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
# Prepare imagery
##################################################################

defol_col = (ee.ImageCollection(f'projects/{args.project}/assets/defoliation_score_North_America')
             .filter(ee.Filter.eq('year', 2021)));
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

forest_change = ee.Image("UMD/hansen/global_forest_change_2024_v1_12")
forest_change_mask = forest_change.select('lossyear').lte(2021).unmask().Not()


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

    if args.submit:
        task = ee.batch.Export.image.toDrive(
            image=defol_tile,
            description=f'upscaled_{i}',
            folder='Upscaled_2021',
            region=gridCell, 
            scale=3000,
            crs='EPSG:5070',
            maxPixels=1e10
        )
        task.start()