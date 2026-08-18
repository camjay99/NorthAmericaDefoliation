import argparse
import json

import ee

import geometries
import preprocessing
import submission


##############################################################
# Parse arguments
##############################################################

parser = argparse.ArgumentParser(
    description='Options for calculating maximum early-season EVI across all years')

# The script will ONLY submit the run when -s or --submit is included.
parser.add_argument('--submit', '-s', action='store_true')

# The script will ONLY create image manifests when -m or --create_manifests is included.
parser.add_argument('--create_manifests', '-i', action='store_true')

# Whether to export results to a cloud storage bucket. If true,
# `bucket` must also be set.
parser.add_argument('--cloudstorage', '-C', action='store_true')

# Cloud storage bucket to save results in.
parser.add_argument('--bucket', '-b', action='store', default=None)

# The project to submit the code in.
# You may be prompted to to authenticate.
parser.add_argument('--project', '-p', action='store',
                    default=None, required=True)

# The first and last years over which to compute the maximum early-season EVI.
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
    file_name_prefix = f'max_early_season_evi_{name}/max_early_season_evi_all_years_v1_{args.data}'
    image_manifests = {}
assetID = f'projects/{args.project}/assets/max_early_season_evi_{name}/max_early_season_evi_all_years_v1_{args.data}'
description_base = f'{name}_MaxEarlySeasonEVI_AllYears_{args.data}'

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

for i in range(gridSize):
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

    # Pull imagery across the full range of years at once so the maximum
    # early-season EVI is computed over all years rather than per year.
    start_date = ee.Date.fromYMD(args.model_start, 1, 1)
    end_date = ee.Date.fromYMD(args.model_end + 1, 1, 1)

    if args.data == 'HLS':
        col = preprocessing.preprocess_HLS(start_date, end_date,
                                           gridCell, 90, 200,
                                           phenology=early_season)


    ##########################################
    # Calculate maximum observed early-season EVI across all years
    ##########################################

    max_evi = col.select('EVI').max().rename('max_EVI')
    max_evi = (max_evi.set('source', args.data)
                    .set('model_start', args.model_start)
                    .set('model_end', args.model_end)
                    .set('window', args.window)
                    .set('min', args.min)
                    .set('max', args.max)
                    .set('project', 'NorthAmerica'))
    max_evi = (max_evi.subtract(args.min)
                    .divide(args.max-args.min).multiply(65_535).uint16())


    #################################
    # Submit batch job
    #################################

    if args.submit:
        submission.submit_job(
            image=max_evi,
            assetID=assetID,
            file_name_prefix=file_name_prefix,
            description_base=description_base,
            scale=preprocessing.resolutions[args.data],
            crs=args.crs,
            region=gridCell,
            cloudstorage=args.cloudstorage,
            bucket=args.bucket,
            i=i
        )
    if args.create_manifests:
        image_manifests[i] = submission.create_manifest(
            assetID=assetID,
            file_name_prefix=file_name_prefix,
            description_base=description_base,
            model_start=args.model_start,
            model_end=args.model_end,
            properties={
                'source': args.data,
                'model_start': args.model_start,
                'model_end': args.model_end,
                'window': args.window,
                'min': args.min,
                'max': args.max,
                'project': 'NorthAmerica'
            },
            bucket=args.bucket,
            i=i
        )
if args.cloudstorage:
    with open("image_manifests.json", 'w')  as f:
        json.dump(image_manifests, f)
