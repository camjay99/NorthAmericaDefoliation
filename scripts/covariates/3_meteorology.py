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

# The first and last years to look for calculating monthly meteorology.
parser.add_argument('--start', '-S', action='store', type=int, default=2019)
parser.add_argument('--end', '-E', action='store', type=int, default=2023)

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
file_name_prefix = f'meteorology_North_America/meteorology'
image_manifests = {}
description_base = f'North_America_Meteorology'
geometry = geometries.get_range()


##################################################################
# Iterate over each year/month and calculate average climate
##################################################################

for year in range(args.start, args.end + 1):
    for i, month in enumerate(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
                               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']):
        ##################################################################
        # Prepare Data
        ##################################################################
        
        met = (ee.ImageCollection('NASA/ORNL/DAYMET_V4')
                 .select(["tmin", "tmax", "prcp"])
                 .filterDate(f'{year}-01-01', f'{year+1}-01-01')
                 .filter(ee.Filter.calendarRange(i + 1, i + 1, 'month')))

        proj = met.first().projection().getInfo()


        temp = (met
                  .select(['tmin', 'tmax'])
                  .mean()
                  .add(60)
                  .multiply(65535./120)
                  .toUint16())
        prcp = (met
                  .select(['prcp'])
                  .sum()
                  .multiply(65535./16864)
                  .toUint16())
        met = ee.Image([temp, prcp])

        ##################################################
        # Export results
        ##################################################

        if args.submit:
            submission.submit_job(
                image=met,
                assetID='',
                file_name_prefix=f'{file_name_prefix}_{year}_{month}',
                description_base=f'{description_base}_{year}_{month}',
                crsTransform=proj['transform'],
                crs=proj['wkt'],
                region=geometry,
                cloudstorage=True,
                bucket=args.bucket,
            )