import os

import requests

DATA_API_URL = 'https://api.planet.com/data/v1'
QUICK_SEARCH_URL = f'{DATA_API_URL}/quick-search'

ORDERS_API_URL = 'https://api.planet.com/compute/ops/orders/v2'

# Default mapping from Data API asset types to Orders API product bundles.
# Override with --product-bundle if a scene list uses a different asset type.
ASSET_TYPE_TO_PRODUCT_BUNDLE = {
    'ortho_analytic_4b': 'analytic_udm2',
    'ortho_analytic_4b_sr': 'analytic_sr_udm2',
    'ortho_analytic_8b': 'analytic_8b_udm2',
    'ortho_analytic_8b_sr': 'analytic_8b_sr_udm2',
    'ortho_visual': 'visual',
}

def get_api_key(api_key=None):
    api_key = api_key or os.environ.get('PL_API_KEY')
    assert api_key, (
        'Planet API key not provided. Pass api_key or set the '
        'PL_API_KEY environment variable.')
    return api_key


def bbox_to_geometry(bbox):
    # bbox is [min_lon, min_lat, max_lon, max_lat]
    min_lon, min_lat, max_lon, max_lat = bbox
    return {
        'type': 'Polygon',
        'coordinates': [[
            [min_lon, min_lat],
            [max_lon, min_lat],
            [max_lon, max_lat],
            [min_lon, max_lat],
            [min_lon, min_lat],
        ]],
    }


def build_filter(geometry, start_date=None, end_date=None,
                  min_cloud_cover=None, max_cloud_cover=None,
                  asset_type='ortho_analytic_4b'):
    filters = [{
        'type': 'GeometryFilter',
        'field_name': 'geometry',
        'config': geometry,
    }]

    if start_date or end_date:
        date_range = {}
        if start_date:
            date_range['gte'] = start_date
        if end_date:
            date_range['lte'] = end_date
        filters.append({
            'type': 'DateRangeFilter',
            'field_name': 'acquired',
            'config': date_range,
        })

    if min_cloud_cover is not None or max_cloud_cover is not None:
        cloud_range = {}
        if min_cloud_cover is not None:
            cloud_range['gte'] = min_cloud_cover
        if max_cloud_cover is not None:
            cloud_range['lte'] = max_cloud_cover
        filters.append({
            'type': 'RangeFilter',
            'field_name': 'cloud_cover',
            'config': cloud_range,
        })

    if asset_type:
        filters.append({
            'type': 'AssetFilter',
            'config': [asset_type],
        })

    return {'type': 'AndFilter', 'config': filters}


def search(api_key, item_types, search_filter, limit=None):
    session = requests.Session()
    session.auth = (api_key, '')

    body = {'item_types': item_types, 'filter': search_filter}
    url = QUICK_SEARCH_URL

    results = []
    while url:
        if url == QUICK_SEARCH_URL:
            response = session.post(url, json=body)
        else:
            # Pagination links already carry the search parameters.
            response = session.get(url)
        response.raise_for_status()
        page = response.json()

        results.extend(page['features'])
        if limit is not None and len(results) >= limit:
            return results[:limit]

        url = page.get('_links', {}).get('_next')

    return results


def scenes_to_records(scenes, asset_type='ortho_analytic_4b'):
    return [{
        'id': scene['id'],
        'item_type': scene['properties']['item_type'],
        'asset_type': asset_type,
        'acquired': scene['properties'].get('acquired'),
        'cloud_cover': scene['properties'].get('cloud_cover'),
        'geometry': scene['geometry'],
    } for scene in scenes]


def product_bundle_for_asset_type(asset_type):
    product_bundle = ASSET_TYPE_TO_PRODUCT_BUNDLE.get(asset_type)
    assert product_bundle, (
        f'No known product bundle for asset type {asset_type!r}. '
        'Pass --product-bundle explicitly.')
    return product_bundle


def build_order(name, item_type, item_ids, product_bundle, archive_type=None):
    product = {
        'item_ids': item_ids,
        'item_type': item_type,
        'product_bundle': product_bundle,
    }

    order = {'name': name, 'products': [product]}

    if archive_type:
        order['delivery'] = {
            'archive_type': archive_type,
            'archive_filename': '{{name}}_{{order_id}}.zip',
        }

    return order


def submit_order(api_key, order):
    session = requests.Session()
    session.auth = (api_key, '')

    response = session.post(ORDERS_API_URL, json=order)
    response.raise_for_status()
    return response.json()
