import argparse
import json

import ee

import geometries
import preprocessing


##############################################################
# Parse arguments
##############################################################

parser = argparse.ArgumentParser(
    description='Options for calculating maximum early-season EVI')

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

# The first and last years to calculate the maximum early-season EVI for.
parser.add_argument('--start', '-S', action='store', type=int, default=2019)
parser.add_argument('--end', '-E', action='store', type=int, default=2023)

# The first and last years of the average phenology model used to
# identify the start of the growing season (SoS).
parser.add_argument('--model_start', '-n', action='store', type=int, default=2019)
parser.add_argument('--model_end', '-N', action='store', type=int, default=2023)

# Number of days after SoS to include when searching for the maximum EVI.
parser.add_argument('--window', '-e', action='store', type=int, default=30)

# The data source to use.
parser.add_argument('--data', '-d', action='store',
                    default='HLS', choices=preprocessing.sources)

# The geomtry to calculate the maximum early-season EVI within.
# A list of valid geometries are available in scripts/geometries.py
parser.add_argument('--geometry', '-g', action='store',
                    default=None, choices=geometries.site_names)

# State to calculate the maximum early-season EVI within.
parser.add_argument('--state', '-x', action='store',
                    default=None)

# Use total Spongy Moth Range to calculate the maximum early-season EVI within.
parser.add_argument('--range', '-R', action='store_true')

# The CRS to output the resulting layers in.
parser.add_argument('--crs', '-c', action='store', default='epsg:5070')

# The width/length of grid cells to use for computation (in lat/lon degrees)
parser.add_argument('--width', '-w', action='store', type=float, default=0.75)
parser.add_argument('--length', '-l', action='store', type=float, default=0.75)

# The min/max of EVI used for data compression
parser.add_argument('--min', '-m', action='store', type=float, default=0)
parser.add_argument('--max', '-M', action='store', type=float, default=1)

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
    file_name_prefix = f'max_early_season_evi_{name}/max_early_season_evi_v1_{args.data}'
    image_manifests = {}
assetID = f'projects/{args.project}/assets/max_early_season_evi_{name}/max_early_season_evi_v1_{args.data}'
description_base = f'{name}_MaxEarlySeasonEVI_{args.data}'

pheno_coll = ee.ImageCollection(f'projects/{args.project}/assets/average_phenology_{name}')
pheno_coll = (pheno_coll.filter(ee.Filter.eq('source', args.data))
                        .filter(ee.Filter.eq('start', args.model_start))
                        .filter(ee.Filter.eq('end', args.model_end)))


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

for i in range(331, gridSize):
    gridCell = ee.Feature(gridList.get(i)).geometry()

    ##################################################################
    # Prepare Data
    ##################################################################

    # Use maximum composite as there are some oddities in overlapping
    # gridCells with some corners missing data.
    phenology = pheno_coll.filterBounds(gridCell).max()

    # Restrict the phenology window to the start of the growing season:
    # from SoS through SoS + window days.
    early_season = ee.Image([
        phenology.select('SoS'),
        phenology.select('SoS').add(args.window).rename('EoS')
    ])

    years = list(range(args.start, args.end + 1))

    for year in years:
        start_date = ee.Date.fromYMD(year, 1, 1)
        end_date = ee.Date.fromYMD(year + 1, 1, 1)

        if args.data == 'HLS':
            col = preprocessing.preprocess_HLS(start_date, end_date,
                                               gridCell, 90, 200,
                                               phenology=early_season)


        ##########################################
        # Calculate maximum observed early-season EVI
        ##########################################

        max_evi = col.select('EVI').max().rename('max_EVI')
        max_evi = (max_evi.set('source', args.data)
                        .set('model_start', args.model_start)
                        .set('model_end', args.model_end)
                        .set('window', args.window)
                        .set('min', args.min)
                        .set('max', args.max)
                        .set('year', year)
                        .set('project', 'NorthAmerica'))
        max_evi = (max_evi.subtract(args.min)
                        .divide(args.max-args.min).multiply(65_535).uint16())


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
                    image=max_evi,
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
                        'model_start':args.model_start,
                        'model_end':args.model_end,
                        'window':args.window,
                        'max':args.max,
                        'min':args.min,
                        'year':year,
                        'project':'NorthAmerica'
                    },
                    'tilesets': [
                        {'id': '0', 'sources': [ {'uris': [f'gs://{args.bucket}/{image_name}.tif']}]}
                    ],
                    'startTime': f'{year}-01-01T00:00:00.000000000Z',
                    'endTime': f'{year+1}-01-01T00:00:00.000000000Z'
                }
            else:
                imageName = f'{assetID}_{year}'
                description = f'{description_base}_{year}'
                if gridSize > 1:
                    imageName = f'{assetID}_{year}_tile_{i}'
                    description = f'{description}_tile_{i}'

                task = ee.batch.Export.image.toAsset(
                    image=max_evi,
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
