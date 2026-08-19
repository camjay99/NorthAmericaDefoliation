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
    description='Options for calculating late-season defoliation')

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

# The min/max of defoliation for data compression
parser.add_argument('--min', '-m', action='store', type=float, default=-1)
parser.add_argument('--max', '-M', action='store', type=float, default=1)

# Day of year after which to start looking for late-season declines. 
parser.add_argument('--start_doy', '-D', action='store', type=int, default=210)

# Hard upper bound on the day of year used to filter the input collection,
# wide enough to include any real end of season (EoS). The phenology
# model's per-pixel EoS still clips the true window.
parser.add_argument('--end_doy', '-e', action='store', type=int, default=365)

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
    file_name_prefix = f'late_season_defoliation_score_{name}/late_season_defoliation_score_v1_{args.data}'
    image_manifests = {}
assetID = f'projects/{args.project}/assets/late_season_defoliation_score_{name}/late_season_defoliation_score_v1_{args.data}'
description_base = f'{name}_LateSeasonDefoliation_{args.data}'

pheno_coll = ee.ImageCollection(f'projects/{args.project}/assets/average_phenology_{name}')
pheno_coll = (pheno_coll.filter(ee.Filter.eq('source', args.data))
                        .filter(ee.Filter.eq('start', args.model_start))
                        .filter(ee.Filter.eq('end', args.model_end)))
model_coll = (ee.ImageCollection(f'projects/{args.project}/assets/seasonal_trend_{name}')
              .filter(ee.Filter.eq('start', args.model_start))
              .filter(ee.Filter.eq('end', args.model_end))
              .filter(ee.Filter.eq('method', 'Theil-Sen'))
              .filter(ee.Filter.eq('source', args.data))
              .filter(ee.Filter.eq('project', 'NorthAmerica')))
swir_coll = (ee.ImageCollection(f'projects/{args.project}/assets/swir_summary_{name}')
             #.filter(ee.Filter.eq('start', args.model_start))
             #.filter(ee.Filter.eq('end', args.model_end))
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
    gridCell = ee.Feature(gridList.get(i)).geometry()

    ##################################################################
    # Prepare Data
    ##################################################################

    # Use maximum composite as there are some oddities in overlapping
    # gridCells with some corners missing data.
    phenology = pheno_coll.filterBounds(gridCell).max()
    models = model_coll.filterBounds(gridCell)
    swir_summary = swir_coll.filterBounds(gridCell)

    # Decompress models
    def decompress_model(image):
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
    models = models.map(decompress_model).mosaic().unmask()

    # Decompress SWIR summary
    def decompress_swir(image):
        image = image.updateMask(image) # Mask 0s, as sometimes the data is available in overlapping tiles
        min_swir = ee.Number(image.get('min'))
        max_swir = ee.Number(image.get('max'))
        image = (image.float()
                      .multiply(max_swir.subtract(min_swir).divide(65_535))
                      .add(min_swir))
        return image
    swir_summary = swir_summary.map(decompress_swir).mosaic().unmask()

    # Compute defoliation score
    years = list(range(args.start, args.end + 1))

    for year in years:
        start_date = ee.Date.fromYMD(year, 1, 1)
        end_date = ee.Date.fromYMD(year + 1, 1, 1)

        if args.data == 'HLS':
            # phenology carries this pixel's real SoS/EoS bands, so the
            # per-pixel mask in preprocess_HLS clips the window to
            # [start_doy, EoS] even though end_doy is a wide fixed bound.
            col = preprocessing.preprocess_HLS(start_date, end_date,
                                               gridCell,
                                               args.start_doy, args.end_doy,
                                               phenology)


        ######################################
        # Estimate defoliation in given window
        ######################################

        # Remove SWIR1 outliers
        def mask_swir_outliers(image):
            # Must be an outlier in both SWIR1 and SWIR2.
            swir1 = image.select('SWIR1')
            swir1_mask = swir1.lte(swir_summary.select('mean_SWIR1')
                            .subtract(swir_summary.select('stddev_SWIR1').multiply(3)))
            swir2 = image.select('SWIR2')
            swir2_mask = swir2.lte(swir_summary.select('mean_SWIR2')
                            .subtract(swir_summary.select('stddev_SWIR2').multiply(3)))
            return image.updateMask(swir1_mask.And(swir2_mask).Not())


        col = col.map(mask_swir_outliers)

        # Calculate anomaly
        def calc_anom(image):
            slope = models.select('slope')
            offset = models.select('offset')
            doy = image.select('doy')
            predict = slope.multiply(doy).add(offset)
            anom = image.select('EVI').subtract(predict)

            return image.addBands(anom.rename('EVI_anom'))

        def calc_statistics(images):
            images = images.map(calc_anom)
            mean_intensity = images.select("EVI_anom").mean().rename("mean_intensity")

            return mean_intensity

        defol = calc_statistics(col)
        defol = (defol.set('source', args.data)
                      .set('start', args.model_start)
                      .set('end', args.model_end)
                      .set('min', args.min)
                      .set('max', args.max)
                      .set('year', year)
                      .set('start_doy', args.start_doy)
                      .set('end_doy', args.end_doy)
                      .set('project', 'NorthAmerica')
                      .set('method', 'Theil-Sen'))
        defol = (defol.subtract(args.min)
                    .divide(args.max-args.min).multiply(65_535).uint16())


        #################################
        # Submit batch job
        #################################

        if args.submit:
            submission.submit_job(
                image=defol,
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
                    'max':args.max,
                    'min':args.min,
                    'year':year,
                    'start_doy':args.start_doy,
                    'end_doy':args.end_doy,
                    'project':'NorthAmerica',
                    'method':'Theil-Sen'
                },
                bucket=args.bucket,
                i=i
            )
if args.cloudstorage:
    with open("image_manifests.json", 'w')  as f:
        json.dump(image_manifests, f)
