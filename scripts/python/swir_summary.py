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
parser.add_argument('--start', '-S', action='store', type=int, default=2019)
parser.add_argument('--end', '-E', action='store', type=int, default=2023)

# The data source to use for calculating trends.
parser.add_argument('--data', '-d', action='store', 
                    default='HLS', choices=preprocessing.sources)

# The geomtry to calculate defoliation within. 
# A list of valid geometries are available in scripts/geometries.py
parser.add_argument('--geometry', '-g', action='store', 
                    default=None, choices=geometries.site_names)

# State to calculate trends within.
parser.add_argument('--state', '-x', action='store', 
                    default=None)

# Use total Spongy Moth Range to calculate trends within.
parser.add_argument('--range', '-R', action='store_true')

# The CRS to output the resulting layers in.
parser.add_argument('--crs', '-c', action='store', default='epsg:5070')

# The width/length of grid cells to use for computation (in lat/lon degrees)
parser.add_argument('--width', '-w', action='store', type=float, default=0.75)
parser.add_argument('--length', '-l', action='store', type=float, default=0.75)

# The min/max of defoliation for data compression
parser.add_argument('--min', '-m', action='store', type=float, default=0)
parser.add_argument('--max', '-M', action='store', type=float, default=0.5)

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

num_not_specified = (args.geometry is None) + (args.state is None) + (not args.range)
assert num_not_specified != 3, "Must specify at least one of geometry, state, or range."
assert num_not_specified > 1, "Only one of geometry, state, and range can be specified at a time."

if args.geometry:
    name = args.geometry
    geometry = geometries.get_geometry(args.geometry)
elif args.state:
    name = args.state.replace(" ", "_")
    geometry = geometries.get_state(args.state)
else: 
    name = "North_America"
    geometry = geometries.get_range()

if args.cloudstorage:
    assert (args.bucket is not None), "Must specify bucket if exporting to cloud storage."
    file_name_prefix = f'swir_summary_{name}/swir_summary_{args.data}'
    image_manifests = {}
assetID = f'projects/{args.project}/assets/swir_summary_{name}/swir_summary_{args.data}'
description_base = f'{name}_SWIR_Summary_{args.data}'


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
    gridCell = ee.Feature(gridList.get(i)).geometry()

    ##################################################################
    # Prepare Data
    ##################################################################

    # Identify SWIR outliers based on all year average.
    if args.data == 'HLS':
        start_date = ee.Date.fromYMD(args.start, 1, 1)
        end_date = ee.Date.fromYMD(args.end + 1, 1, 1)
        full_col = preprocessing.preprocess_HLS(start_date, end_date,
                                                gridCell, 
                                                start_doy=160, end_doy=210)
    
    mean_swir = full_col.select(['SWIR1', 'SWIR2']).mean()
    std_dev_swir = full_col.select(['SWIR1', 'SWIR2']).reduce(ee.Reducer.stdDev())


    swir_summary = ee.Image([mean_swir.rename(['mean_SWIR1', 'mean_SWIR2']),
                            std_dev_swir.rename(['stddev_SWIR1', 'stddev_SWIR2'])])
    swir_summary = (swir_summary.set('source', args.data)
                    .set('min', args.min)
                    .set('max', args.max)
                    .set('project', 'NorthAmerica')
                    .set('method', 'SWIR Summary'))
    swir_summary = (swir_summary.subtract(args.min)
                .divide(args.max-args.min).multiply(65_535).uint16())


    #################################
    # Submit batch job
    #################################

    if args.submit:
        if args.cloudstorage:
            # Save in a Cloud Storage Bucket
            if gridSize > 1:
                asset_name = f'{assetID}_tile_{i}'
                image_name = f'{file_name_prefix}_tile_{i}'
                description = f'{description_base}_tile_{i}'
            else:
                asset_name = assetID
                image_name = file_name_prefix
                description = description_base

            # task = ee.batch.Export.image.toCloudStorage(
            #     image=swir_summary,
            #     description=description,
            #     bucket=args.bucket,
            #     fileNamePrefix=image_name,
            #     region=gridCell,
            #     scale=preprocessing.resolutions[args.data],
            #     crs=args.crs,
            #     maxPixels=1e10,
            #     formatOptions={
            #         'cloudOptimized': True,
            #     }
            # )
            # task.start()
            
            # Create an image manifest for adding image as an asset
            image_manifests[f"{i}"] = {
                'name': asset_name,
                'properties': {
                    'source':args.data,
                    'max':args.max,
                    'min':args.min,
                    'project':'NorthAmerica',
                    'method':'SWIR Summary'
                },
                'tilesets': [
                    {'id': '0', 'sources': [ {'uris': [f'gs://{args.bucket}/{image_name}.tif']}]}
                ],
                'startTime': f'{args.start}-01-01T00:00:00.000000000Z',
                'endTime': f'{args.end+1}-01-01T00:00:00.000000000Z'
            }
        else:
            imageName = assetID
            description = description_base
            if gridSize > 1:
                imageName = f'{assetID}_tile_{i}'
                description = f'{description}_tile_{i}'

            task = ee.batch.Export.image.toAsset(
                image=swir_summary,
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