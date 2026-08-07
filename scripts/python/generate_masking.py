import argparse
import json

import ee

import geometries
import preprocessing

##############################################################
# Parse arguments
##############################################################


parser = argparse.ArgumentParser(
    description='Options for calculating seasonal trends')

# The script will ONLY submit the run when -s or --submit is included.
parser.add_argument('--submit', '-s', action='store_true')

# Whether to export results to a cloud storage bucket. If true,
# `bucket` must also be set.
parser.add_argument('--cloudstorage', '-C', action='store_true')

# Cloud storage bucket to save results in.
parser.add_argument('--bucket', '-b', action='store', default=None)

# The project to submit the code in. You may be prompted to to authenticate.
parser.add_argument('--project', '-p', action='store', default=None, required=True)

# The first and last years to look for defoliation signals in.
parser.add_argument('--start', '-S', action='store', type=int, default=2019)
parser.add_argument('--end', '-E', action='store', type=int, default=2023)

# The first and last years of baseline models used.
parser.add_argument('--model_start', '-n', action='store', type=int, default=2019)
parser.add_argument('--model_end', '-N', action='store', type=int, default=2023)

# The data source to use for calculating trends.
parser.add_argument('--data', '-d', action='store', 
                    default='HLS', choices=preprocessing.sources)

# The geomtry to calculate defoliation within. A list of valid geometries are available in scripts/geometries.py
parser.add_argument('--geometry', '-g', action='store', default='Mt_Pleasant', choices=geometries.site_names)

# State to calculate defoliation over. If specified, geometry is ignored
parser.add_argument('--state', '-t', action='store', default=None)

# Use total Spongy Moth Range to calculate trends within.
parser.add_argument('--range', '-R', action='store_true')

# The geomtry to calculate defoliation within. A list of valid geometries are available in scripts/geometries.py
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
# Specify base names and load previous results for small site
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
    file_name_prefix = f'qa_mask_{name}/qa_mask_{args.data}'
    image_manifests = {}
    if args.rescale:
        file_name_prefix += '_rescaled'
assetID = f'projects/{args.project}/assets/qa_mask_{name}/qa_mask_{args.data}'
description_base = f'{name}_qa_mask_{args.data}'


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

    years = list(range(args.start, args.end + 1))
    for year in years:
        year_start_date = ee.Date.fromYMD(year, 1, 1)
        year_end_date = ee.Date.fromYMD(year + 1, 1, 1)
        start_date = ee.Date.fromYMD(args.model_start, 1, 1)
        end_date = ee.Date.fromYMD(args.model_end, 1, 1)
        # Collect observations
        if args.data == 'HLS':
            year_col = preprocessing.preprocess_HLS(year_start_date, year_end_date,
                                                    gridCell, None, False)
            all_col = preprocessing.preprocess_HLS(start_date, end_date,
                                                gridCell, None, False)
            
        # Create yearly observation mask
        def yearly_obs_mask(col, threshold):
            col = col.filter(ee.Filter.dayOfYear(161, 208))

            # Combine days with multipe observations
            withDates = col.map(lambda image: image.set('date', image.date().format('YYYY-MM-dd')))

            mosaicList = (withDates.aggregate_array('date')
                .distinct()
                .map(lambda date: col.filterDate(ee.Date(date), ee.Date(date).advance(1, 'day')).max()))
                
            obs_counts = (ee.ImageCollection.fromImages(mosaicList)
                .sum()
                .gte(threshold)
                .toUint16()
                .unmask(0))
            
            return obs_counts

        strong_obs_mask = yearly_obs_mask(year_col, 3)

        weak_obs_mask = yearly_obs_mask(year_col, 2)

        # Preseason max across all years
        preseason_max_all_years = (all_col.filter(ee.Filter.dayOfYear(130, 170))
                                    .select('EVI')
                                    .reduce(ee.Reducer.percentile([95])))

        ## Create yearly preseason mask
        def yearly_preseason_mask(col, threshold):  
            col = col.filter(ee.Filter.dayOfYear(130, 170)).select('EVI')
            preseason_max_year = col.reduce(ee.Reducer.percentile([95]))
            preseason_max_count = (col.reduce(ee.Reducer.count())
                                        .eq(0)
                                        .unmask(0))
            preseason_gap = (preseason_max_year.subtract(preseason_max_all_years)
                                                .gte(threshold)
                                                .unmask(0))

            # Ensure we have large gap and observtions to base this on.
            return preseason_gap.Or(preseason_max_count).toUint16()

        preseason_mask = yearly_preseason_mask(year_col, -0.15)

        # Create forest cover mask
        ## Load NLCD 2019 landcover map
        nlcd_landcover = ee.ImageCollection('USGS/NLCD_RELEASES/2019_REL/NLCD') \
            .filter(ee.Filter.eq('system:index', '2019')).first().select('landcover')
        forest_mask = nlcd_landcover.gte(41).And(nlcd_landcover.lte(43)).toUint16()

        # Combine masks into single image
        qa_mask = (ee.Image([forest_mask.rename('forest'),
                            preseason_mask.rename('preseason'),
                            strong_obs_mask.rename('count_3'),
                            weak_obs_mask.rename('count_2')])
                    .set('source', args.data)
                    .set('start', args.model_start)
                    .set('end', args.model_end)
                    .set('year', year)
                    .set('project', 'NorthAmerica'))
        
        # Create Task
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
                    image=qa_mask,
                    description=description,
                    bucket=args.bucket,
                    fileNamePrefix=image_name,
                    region=gridCell,
                    scale=preprocessing.resolutions[args.data],
                    crs=args.crs,
                    maxPixels=1e10,
                    formatOptions={
                        'cloudOptimized': True,
                    },
                    pyramidingPolicy = {'.default': 'sample'},
                )
                task.start()

                # Create an image manifest for adding image as an asset
                image_manifests[f"{year}_{i}"] = {
                    'name': asset_name,
                    'properties': {
                        'source':args.data,
                        'start':args.model_start,
                        'end':args.model_end,
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
                if gridSize > 1:
                    image_name = f'{assetID}_{year}_tile_{i}'
                    description = f'{description_base}_{year}_tile_{i}'
                else:
                    image_name = f'{assetID}_{year}'
                    description = f'{description_base}_{year}'
                task = ee.batch.Export.image.toAsset(
                    image=qa_mask,
                    description=description,
                    assetId=image_name,
                    region=gridCell, 
                    scale=preprocessing.resolutions[args.data],
                    crs=args.crs,
                    pyramidingPolicy={'.default': 'sample'},
                    maxPixels=1e10
                )
                task.start()
if args.cloudstorage:
    with open("image_manifests.json", 'w')  as f:
        json.dump(image_manifests, f)