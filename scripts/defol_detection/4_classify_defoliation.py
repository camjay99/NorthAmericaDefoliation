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
    description='Options for classifying defoliation status')

# The script will ONLY submit the run when -s or --submit is included.
parser.add_argument('--submit', '-s', action='store_true')

# The script will ONLY create image manifests when -i or --create_manifests is included.
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

# The first and last years to look for defoliation signals in.
parser.add_argument('--start', '-S', action='store', type=int, default=2019)
parser.add_argument('--end', '-E', action='store', type=int, default=2023)

# The first and last years of baseline models used.
parser.add_argument('--model_start', '-n', action='store', type=int, default=2019)
parser.add_argument('--model_end', '-N', action='store', type=int, default=2023)

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
    file_name_prefix = f'forest_status_{name}/forest_status_v1_{args.data}'
    image_manifests = {}
assetID = f'projects/{args.project}/assets/forest_status_{name}/forest_status_v1_{args.data}'
description_base = f'{name}_ForestStatus_{args.data}'

all_years_early_season_coll = (
    ee.ImageCollection(
        f'projects/{args.project}/all_years_max_early_season_evi_{name}/')
        .filter(ee.Filter.eq('source', args.data))
        .filter(ee.Filter.eq('start', args.model_start))
        .filter(ee.Filter.eq('end', args.model_end))
        .filter(ee.Filter.eq('project', 'NorthAmerica')))
yearly_early_season_coll = (
    ee.ImageCollection(
        f'projects/{args.project}/yearly_max_early_season_evi_{name}/')
        .filter(ee.Filter.eq('source', args.data))
        .filter(ee.Filter.eq('start', args.model_start))
        .filter(ee.Filter.eq('end', args.model_end))
        .filter(ee.Filter.eq('project', 'NorthAmerica')))
mid_season_coll = (
    ee.ImageCollection(
        f'projects/{args.project}/defoliation_score_{name}/')
        .filter(ee.Filter.eq('source', args.data))
        .filter(ee.Filter.eq('start', args.model_start))
        .filter(ee.Filter.eq('end', args.model_end))
        .filter(ee.Filter.eq('project', 'NorthAmerica'))
        .filter(ee.Filter.eq('method', 'Theil-Sen')))
late_season_coll = (
    ee.ImageCollection(
        f'projects/{args.project}/assets/late_season_defoliation_score_{name}/')
        .filter(ee.Filter.eq('source', args.data))
        .filter(ee.Filter.eq('start', args.model_start))
        .filter(ee.Filter.eq('end', args.model_end))
        .filter(ee.Filter.eq('project', 'NorthAmerica'))
        .filter(ee.Filter.eq('method', 'Theil-Sen')))


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

    for year in range(args.start, args.end + 1):
        def decompress(image):
            image = image.updateMask(image)
            min = ee.Number(image.get('min'))
            max = ee.Number(image.get('max'))
            image = (image.float()
                          .multiply(max.subtract(min.divide(65_535)))
                          .add(min))
            return image
        # Use maximum composite as there are some oddities in overlapping
        # gridCells with some corners missing data.
        all_years_es = (all_years_early_season_coll
                            .filterBounds(gridCell)
                            .map(decompress)
                            .mosaic()
                            .unmask())
        yearly_es = (yearly_early_season_coll
                     .filterBounds(gridCell)
                     .filter(ee.Filter.eq('year', year))
                     .map(decompress)
                     .mosaic()
                     .unmask())
        mid_season = (mid_season_coll
                      .filterBounds(gridCell)
                      .filter(ee.Filter.eq('year', year))
                      .map(decompress)
                      .mosaic()
                      .unmask())
        late_season = (late_season_coll
                       .filterBounds(gridCell)
                       .filter(ee.Filter.eq('year', year))
                       .map(decompress)
                       .mosaic()
                       .unmask())

        # Create classifications
        es_dec = yearly_es.subtract(all_years_es).lte(-0.3)
        ms_dec = mid_season.lte(-0.04)
        ls_dec = late_season.lte(-0.04)

        status = (es_dec.multiply(4)
                   .add(ms_dec.multiply(2))
                   .add(ls_dec.multiply(1))
                   .toUint8()
                   .rename('forest_status'))
        status = (status.set('source', args.data)
                        .set('start', args.model_start)
                        .set('end', args.model_end)
                        .set('project', 'NorthAmerica')
                        .set('year', year))


        #################################
        # Submit batch job
        #################################

        if args.submit:
            submission.submit_job(
                image=status,
                assetID=assetID,
                file_name_prefix=file_name_prefix,
                description_base=description_base,
                year=year,
                scale=preprocessing.resolutions[args.data],
                crs=args.crs,
                region=gridCell,
                cloudstorage=args.cloudstorage,
                bucket=args.bucket,
                i=i
            )
        if args.create_manifests:
            image_manifests[f"{year}_{i}"] = submission.create_manifest(
                assetID=assetID,
                file_name_prefix=file_name_prefix,
                description_base=description_base,
                year=year,
                properties={
                    'source':args.data,
                    'start':args.model_start,
                    'end':args.model_end,
                    'year':year,
                    'project':'NorthAmerica',
                },
                bucket=args.bucket,
                i=i
            )
if args.cloudstorage:
    with open("image_manifests.json", 'w')  as f:
        json.dump(image_manifests, f)


