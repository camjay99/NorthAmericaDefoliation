import argparse
import json

import ee

import geometries
import masking
import preprocessing


##############################################################
# Parse arguments
##############################################################

parser = argparse.ArgumentParser(
    description='Options for denoising defoliation classifications')

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

# The first and last years to look for defoliation signals in.
parser.add_argument('--start', '-S', action='store', default=2019)
parser.add_argument('--end', '-E', action='store', default=2023)

# The data source to use for calculating trends.
parser.add_argument('--data', '-d', action='store', 
                    default='Sentinel2', choices=preprocessing.sources)

# The geomtry to calculate defoliation within. 
# A list of valid geometries are available in scripts/geometries.py
parser.add_argument('--geometry', '-g', action='store', 
                    default='Mt_Pleasant', choices=geometries.site_names)

# State to calculate trends within.
parser.add_argument('--state', '-x', action='store', 
                    default=None)

# The CRS to output the resulting layers in.
parser.add_argument('--crs', '-c', action='store', default='epsg:5070')

# The width/length of grid cells to use for computation (in lat/lon degrees)
parser.add_argument('--width', '-w', action='store', type=int, default=0.75)
parser.add_argument('--length', '-l', action='store', type=int, default=0.75)

# Rescale within each year to minimize interannual variability
parser.add_argument('--rescale', '-r', action='store_true')

# Select the masking function to use
parser.add_argument('--mask', '-m', action='store', 
                    default='New_York_v1', choices=masking.masks)

# The threshold to use for classification.
parser.add_argument('--threshold', '-t', action='store', type=float, default=-0.040)

# Parse arguments provided to script
args = parser.parse_args()


##############################################################
# Initialize Google Earth Engine API
##############################################################

import ee 

try:
    ee.Initialize(project=args.project)
except:
    # need to authenticate with your credential at the first time
    ee.Authenticate()
    ee.Initialize(project=args.project)


##################################################################
# Specify base names and load previous results
##################################################################

if args.state == None:
    name = args.geometry
    geometry = geometry = geometries.get_geometry(args.geometry)
else:
    name = args.state.replace(" ", "_")
    geometry = geometries.get_state(args.state)

if args.cloudstorage:
    assert (args.bucket is not None), "Must specify bucket if exporting to cloud storage."
    file_name_prefix = f'defoliation_class_{name}/defoliation_class_v2_{args.data}'
    image_manifests = {}
    if args.rescale:
        file_name_prefix += '_rescaled'
assetID = f'projects/{args.project}/assets/defoliation_class_{name}/defoliation_class_v2_{args.data}'
description_base = f'{name}_Classification_{args.data}'
if args.rescale:
    assetID += '_rescaled'
    description_base += '_rescaled'

defol_coll = (ee.ImageCollection(f'projects/{args.project}/assets/defoliation_score_{name}')
              .filter(ee.Filter.eq('rescaled', args.rescale))
              .filter(ee.Filter.eq('start', args.start))
              .filter(ee.Filter.eq('end', args.end))
              .filter(ee.Filter.eq('method', 'Theil-Sen'))
              .filter(ee.Filter.eq('source', args.data))
              .filter(ee.Filter.eq('project', 'NorthAmerica')))


##################################################################
# Split study regions into grid cells of specified size.
##################################################################

#Specify grid size in projection, x and y units (based on projection).
projection = 'EPSG:4326'; # WGS84 lat lon

# Make grid and visualize.
proj = ee.Projection(projection).scale(args.width, args.length)
grid = geometry.coveringGrid(proj)

gridSize = grid.size().getInfo()
gridList = grid.toList(gridSize)

for i in range(gridSize):

    ###############################################################################
    # Create buffer around region to ensure accurate denoising at grid cell edges
    ###############################################################################
    
    gridCell = ee.Feature(gridList.get(i)).geometry()
    gridCellBuffer = gridCell.buffer(distance=400)

    years = list(range(args.start, args.end + 1))

    for year in years:
        ################################################
        # Classify pixels within buffered grid cell
        ################################################

        if args.mask == 'New_York_v1':
            qa_mask = masking.generate_qa_mask_v1(year, gridCellBuffer)
        
        defol_gridCell = (defol_coll.filter(ee.Filter.eq('year', year))
                                    .filterBounds(gridCellBuffer)
                                    .mosaic()
                                    .clip(gridCellBuffer))
        class_gridCell = (defol_gridCell.lte(args.threshold)
                                        .updateMask(qa_mask))
        
        ################################################
        # Remove pixels that are part of small patches
        ################################################

        groups_gridCell = class_gridCell.connectedPixelCount(15, False)
        too_small = groups_gridCell.lte(10)
        
        groups_denoised = class_gridCell.And(too_small.Not()).toUint8()
        
        groups_denoised = (groups_denoised.set('year', year)
                                          .set('rescaled', args.rescale)
                                          .set('start', args.start)
                                          .set('end', args.end)
                                          .set('source', args.data)
                                          .set('project', 'NorthAmerica'))


        #################################
        # Submit batch job
        #################################

        if args.submit:
            if args.cloudstorage:
                # Save in a Cloud Storage Bucket
                if gridSize > 1:
                    asset_name = f'{assetID}_{year}_tile_{i}'
                    image_name = f'{file_name_prefix}_{year}_tile_{i}'
                    description = f'{description_base}_{year}_tile_{i}'
                else:
                    asset_name = f'{assetID}_{year}'
                    image_name = f'{file_name_prefix}_{year}'
                    description = f'{description_base}_{year}'

                task = ee.batch.Export.image.toCloudStorage(
                    image=groups_denoised,
                    description=description,
                    bucket=args.bucket,
                    fileNamePrefix=image_name,
                    region=gridCell,
                    scale=preprocessing.resolutions[args.data],
                    crs=args.crs,
                    maxPixels=1e10,
                    formatOptions={
                        'cloudOptimized': True,
                    }
                )
                task.start()
                
                # Create an image manifest for adding image as an asset
                image_manifests[f"{year}_{i}"] = {
                    'name': asset_name,
                    'properties': {
                        'source':args.data,
                        'start':args.start,
                        'end':args.end,
                        'max':args.max,
                        'min':args.min,
                        'rescaled':str(args.rescale),
                        'year':year,
                        'project':'NorthAmerica'
                    },
                    'tilesets': [
                        {'id': '0', 'sources': [ {'uris': [f'gs://{args.bucket}/{image_name}.tif']}]}
                    ],
                    'startTime': f'{args.start}-01-01T00:00:00.000000000Z',
                    'endTime': f'{args.end+1}-01-01T00:00:00.000000000Z'
                }
            else:
                assetID = f'{assetID}_{year}'
                description = f'{description_base}_{year}'
                if gridSize > 1:
                    imageName = f'{assetID}_tile_{i}'
                    description = f'{description}_tile_{i}'

                task = ee.batch.Export.image.toAsset(
                    image=groups_denoised,
                    description=description,
                    assetId=imageName,
                    region=gridCell, 
                    scale=preprocessing.resolutions[args.data],
                    crs=args.crs,
                    pyramidingPolicy={'.default': 'mean'},
                    maxPixels=1e10
                )
                task.start()
if args.cloudstorage:
    with open("image_manifests.json", 'w')  as f:
        json.dump(image_manifests, f)