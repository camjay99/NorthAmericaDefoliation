import argparse
import json
from collections import defaultdict

import utils.planet as planet

# The Orders API rejects a product with more than this many item_ids, so
# large scene lists are split into multiple orders.
MAX_ITEMS_PER_ORDER = 500


##############################################################
# Parse arguments
##############################################################

parser = argparse.ArgumentParser(
    description='Submit Planet orders for all imagery listed in a scene '
                 'JSON file produced by search_imagery.py.')

parser.add_argument('--input', '-i', action='store', required=True,
                     help='Path to the scene list JSON produced by '
                          'search_imagery.py.')

parser.add_argument('--name', action='store', required=True,
                     help='Base name for the submitted orders. Each order '
                          'is named "{name}_{item_type}_{n}".')

parser.add_argument('--product-bundle', action='store', default=None,
                     help='Orders API product bundle to request for every '
                          'scene. Defaults to a mapping based on each '
                          "scene's asset_type.")

parser.add_argument('--archive-type', action='store', default=None,
                     choices=['zip'],
                     help='If set, deliver each order as a single archive.')

parser.add_argument('--batch-size', action='store', type=int,
                     default=MAX_ITEMS_PER_ORDER,
                     help='Maximum number of scenes per order '
                          f'(max {MAX_ITEMS_PER_ORDER}).')

parser.add_argument('--api-key', action='store', default=None,
                     help='Defaults to the PL_API_KEY environment variable.')

parser.add_argument('--output', '-o', action='store', default=None,
                     help='Path to save the submitted order metadata as '
                          'JSON.')

parser.add_argument('--dry-run', action='store_true',
                     help='Build the orders but do not submit them.')

args = parser.parse_args()

assert 0 < args.batch_size <= MAX_ITEMS_PER_ORDER, (
    f'--batch-size must be between 1 and {MAX_ITEMS_PER_ORDER}.')


##############################################################
# Load scenes and group into orders
##############################################################

with open(args.input) as f:
    records = json.load(f)

assert records, f'No scenes found in {args.input}.'

groups = defaultdict(list)
for record in records:
    if args.product_bundle:
        product_bundle = args.product_bundle
    else:
        product_bundle = planet.product_bundle_for_asset_type(
            record['asset_type'])
    groups[(record['item_type'], product_bundle)].append(record['id'])

print(f'Loaded {len(records)} scenes into {len(groups)} product group(s).')


##############################################################
# Submit orders in batches
##############################################################

api_key = planet.get_api_key(args.api_key)
submitted = []
order_number = 0

for (item_type, product_bundle), item_ids in groups.items():
    for start in range(0, len(item_ids), args.batch_size):
        batch = item_ids[start:start + args.batch_size]
        order_number += 1
        order_name = f'{args.name}_{item_type}_{order_number}'
        order = planet.build_order(
            order_name, item_type, batch, product_bundle,
            archive_type=args.archive_type)

        if args.dry_run:
            print(f'[dry run] {order_name}: {len(batch)} scenes '
                  f'({item_type}, {product_bundle})')
            submitted.append(order)
            continue

        response = planet.submit_order(api_key, order)
        print(f'Submitted {order_name}: {len(batch)} scenes '
              f'({item_type}, {product_bundle}) -> order {response["id"]}')
        submitted.append(response)

print(f'Submitted {len(submitted)} order(s).')


##############################################################
# Save results
##############################################################

if args.output:
    with open(args.output, 'w') as f:
        json.dump(submitted, f, indent=2)
