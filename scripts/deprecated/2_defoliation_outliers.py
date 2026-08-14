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

# Rescale within each year to minimize interannual variability
parser.add_argument('--rescale', '-r', action='store_true')

# The min/max of defoliation for data compression
parser.add_argument('--min', '-m', action='store', type=float, default=-1)
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
    file_name_prefix = f'defoliation_score_{name}/defoliation_score_v3_{args.data}'
    image_manifests = {}
    if args.rescale:
        file_name_prefix += '_rescaled'
assetID = f'projects/{args.project}/assets/defoliation_score_{name}/defoliation_score_v3_{args.data}'
description_base = f'{name}_Defoliation_{args.data}'
if args.rescale:
    assetID += '_rescaled'
    description_base += '_rescaled'

pheno_coll = ee.ImageCollection(f'projects/{args.project}/assets/average_phenology_{name}')
pheno_coll = (pheno_coll.filter(ee.Filter.eq('source', args.data))
                        .filter(ee.Filter.eq('start', args.model_start))
                        .filter(ee.Filter.eq('end', args.model_end)))
model_coll = (ee.ImageCollection(f'projects/{args.project}/assets/seasonal_trend_{name}')
              .filter(ee.Filter.eq('rescaled', str(args.rescale)))
              .filter(ee.Filter.eq('start', args.model_start))
              .filter(ee.Filter.eq('end', args.model_end))
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

for i in [157, 189]:
    gridCell = ee.Feature(gridList.get(i)).geometry()

    ##################################################################
    # Prepare Data
    ##################################################################

    # Use maximum composite as there are some oddities in overlapping 
    # gridCells with some corners missing data.
    phenology = pheno_coll.filterBounds(gridCell).max() 
    models = model_coll.filterBounds(gridCell)
    
    # Uncompress models
    def decompress(image):
        image = image.updateMask(image) # Mask 0s, as sometimes the data is available in overlapping tiles
        min_slope = ee.Number(image.get('min_slope'))
        max_slope = ee.Number(image.get('max_slope'))
        slope = (image.select('slope')
                      .float()
                      .multiply(max_slope.subtract(min_slope).divide(65_535))
                      .add(min_slope))
        min_intercept = ee.Number(image.get('min_intercept'))
        max_intercept = ee.Number(image.get('max_intercept'))
        offset = (image.select('offset')
                       .float()
                       .multiply(max_intercept.subtract(min_intercept).divide(65_535))
                       .add(min_intercept))
        return ee.Image([slope, offset])
    models = models.map(decompress).mosaic().unmask()

    years = list(range(args.start, args.end + 1))

    for year in [2021]:
        start_date = ee.Date.fromYMD(year, 1, 1)
        end_date = ee.Date.fromYMD(year + 1, 1, 1)

        if args.data == 'HLS':
            col = preprocessing.preprocess_HLS(start_date, end_date,
                                               gridCell, phenology)

        if args.rescale:
            col = preprocessing.rescale_years(col, args.start, args.end)

        
        ######################################
        # Estimate defoliation in given window
        ######################################

        # Calculate anomaly
        def calc_anom(image):
            slope = models.select('slope')
            offset = models.select('offset')
            doy = image.select('doy')
            predict = slope.multiply(doy).add(offset)
            anom = image.select('EVI').subtract(predict)

            return image.addBands(anom.rename('EVI_anom'))
        
        def calc_deviation(image, median):
            deviation = image.select('EVI').subtract(median.select('EVI')).abs()
            return deviation.rename('deviation')
        
        def MAD_mask(image, threshold):
            mask = image.select('EVI').gte(threshold)
            return image.updateMask(mask)

        def calc_statistics(images): 
            images = images.filter(ee.Filter.dayOfYear(161, 208))
            images = images.map(calc_anom)
            median = images.median()
            MAD = images.map(lambda image: calc_deviation(image, median)).median().rename('MAD')
            threshold = median.select('EVI').subtract(MAD.multiply(3)).rename('threshold')
            images = images.map(lambda image: MAD_mask(image, threshold))
            mean_intensity = images.select("EVI_anom").mean().rename("mean_intensity")

            return mean_intensity.addBands([median.select('EVI').rename('median'), MAD, threshold])
        
        defol = calc_statistics(col)
        defol = (defol.set('source', args.data)
                      .set('rescaled', args.rescale)
                      .set('start', args.model_start)
                      .set('end', args.model_end)
                      .set('min', args.min)
                      .set('max', args.max)
                      .set('year', year)
                      .set('method', 'Theil-Sen + outliers')
                      .set('project', 'NorthAmerica'))
        defol = (defol.subtract(args.min)
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
                    image=defol,
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
                        'start':args.model_start,
                        'end':args.model_end,
                        'max':args.max,
                        'min':args.min,
                        'rescaled':str(args.rescale),
                        'year':year,
                        'method':'Theil-Sen + outliers',
                        'project':'NorthAmerica'
                    },
                    'tilesets': [
                        {'id': '0', 'sources': [ {'uris': [f'gs://{args.bucket}/{image_name}.tif']}]}
                    ],
                    'startTime': f'{args.start}-01-01T00:00:00.000000000Z',
                    'endTime': f'{args.end+1}-01-01T00:00:00.000000000Z'
                }
            else:
                imageName = f'{assetID}_{year}'
                description = f'{description_base}_{year}'
                if gridSize > 1:
                    imageName = f'{assetID}_tile_{i}'
                    description = f'{description}_tile_{i}'

                task = ee.batch.Export.image.toAsset(
                    image=defol,
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