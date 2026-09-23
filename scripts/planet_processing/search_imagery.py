import argparse
import json

import utils.planet as planet


##############################################################
# Parse arguments
##############################################################

parser = argparse.ArgumentParser(
    description='Search the Planet Data API for 4-band orthorectified '
                 'imagery intersecting a rectangle, and save the results '
                 'for later use placing an order.')

# Rectangle to search, in lon/lat degrees.
parser.add_argument('--bbox', nargs=4, type=float, required=True,
                     metavar=('MIN_LON', 'MIN_LAT', 'MAX_LON', 'MAX_LAT'))

parser.add_argument('--start-date', action='store', default=None,
                     help='ISO 8601 date/time, e.g. 2024-06-01T00:00:00Z')
parser.add_argument('--end-date', action='store', default=None,
                     help='ISO 8601 date/time, e.g. 2024-09-01T00:00:00Z')

parser.add_argument('--min-cloud-cover', action='store', type=float,
                     default=None, help='Fraction between 0 and 1.')
parser.add_argument('--max-cloud-cover', action='store', type=float,
                     default=None, help='Fraction between 0 and 1.')

parser.add_argument('--item-type', action='store', default='PSScene')
parser.add_argument('--asset-type', action='store',
                     default='ortho_analytic_4b')

parser.add_argument('--limit', action='store', type=int, default=None)

parser.add_argument('--api-key', action='store', default=None,
                     help='Defaults to the PL_API_KEY environment variable.')

parser.add_argument('--output', '-o', action='store', required=True,
                     help='Path to save the resulting scene list as JSON.')

args = parser.parse_args()


##############################################################
# Search Planet Data API
##############################################################

api_key = planet.get_api_key(args.api_key)
geometry = planet.bbox_to_geometry(args.bbox)
search_filter = planet.build_filter(
    geometry,
    start_date=args.start_date,
    end_date=args.end_date,
    min_cloud_cover=args.min_cloud_cover,
    max_cloud_cover=args.max_cloud_cover,
    asset_type=args.asset_type,
)

scenes = planet.search(
    api_key, [args.item_type], search_filter, limit=args.limit)
records = planet.scenes_to_records(scenes, asset_type=args.asset_type)

print(f'Found {len(records)} scenes.')


##############################################################
# Save results for later ordering
##############################################################

with open(args.output, 'w') as f:
    json.dump(records, f, indent=2)
