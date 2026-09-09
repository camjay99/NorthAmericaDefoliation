import argparse
import json

import ee
import eeauth

import utils.geometries as geometries
import utils.preprocessing as preprocessing
import utils.submission as submission


##############################################################
# Parse arguments
##############################################################

parser = argparse.ArgumentParser(
    description='Options for calculating growing season')

# The user for submitting jobs.
parser.add_argument('--user', '-u', action='store', default=None)

# The script will ONLY submit the run when -s or --submit is included.
parser.add_argument('--submit', '-s', action='store_true')

# Cloud storage bucket to save results in.
parser.add_argument('--bucket', '-b', action='store', default=None)

# The project to submit the code in. 
# You may be prompted to to authenticate.
parser.add_argument('--project', '-p', action='store', 
                    default=None, required=True)

# Parse arguments provided to script
args = parser.parse_args()


##############################################################
# Initialize Google Earth Engine API
##############################################################

try:
    eeauth.initialize(user=args.user, project=args.project)
except:
    # need to authenticate with your credential at the first time
    eeauth.authenticate(user=args.user)
    eeauth.initialize(user=args.user, project=args.project)


##################################################################
# Specify base names and load previous results
##################################################################

assert (args.bucket is not None), "Must specify bucket if exporting to cloud storage."
file_name_prefix = f'climate_North_America/climate'
image_manifests = {}
description_base = f'North_America_Climate'
geometry = geometries.get_range()


##################################################################
# Iterate over each month and calculate long-term average climate
##################################################################

for i, month in enumerate(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
                           'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
    ##################################################################
    # Prepare Data
    ##################################################################
    
    climate = (ee.ImageCollection('NASA/ORNL/DAYMET_V4')
                 .select(["tmin", "tmax", "prcp"])
                 .filterDate('2000-01-01', '2026-01-01')
                 .filter(ee.Filter.calendarRange(i + 1, i + 1, 'month')))

    proj= climate.first().projection().getInfo()

    temp = (climate
              .select(['tmin', 'tmax'])
              .mean()
              .add(60)
              .multiply(65535./120)
              .toUint16())
    prcp = (climate
              .select(['prcp'])
              .sum()
              .multiply(65535./16864)
              .toUint16())
    climate = ee.Image([temp, prcp])

    ##################################################
    # Export results
    ##################################################

    if args.submit:
        submission.submit_job(
            image=climate,
            assetID='',
            file_name_prefix=f'{file_name_prefix}_{month}',
            description_base=f'{description_base}_{month}',
            crsTransform=proj['transform'],
            crs=proj['wkt'],
            region=geometry,
            cloudstorage=True,
            bucket=args.bucket,

        )