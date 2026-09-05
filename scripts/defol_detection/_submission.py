import time

import ee

def submit_job(
        image,
        assetID,
        file_name_prefix,
        description_base,
        scale,
        crs,
        region,
        cloudstorage,
        bucket=None,
        year=None,
        i=None
):
    asset_name, image_name, description = _create_name_description(
        assetID, file_name_prefix, description_base, year, i
    )

    if cloudstorage:
        # Save in a Cloud Storage Bucket
        task = ee.batch.Export.image.toCloudStorage(
            image=image,
            description=description,
            bucket=bucket,
            fileNamePrefix=image_name,
            region=region,
            scale=scale,
            crs=crs,
            maxPixels=1e10,
            formatOptions={
                'cloudOptimized': True,
            }
        )
        task.start()
    else:
        task = ee.batch.Export.image.toAsset(
            image=image,
            description=description,
            assetId=asset_name,
            region=region,
            scale=scale,
            crs=crs,
            pyramidingPolicy={'.default': 'mean'},
            maxPixels=1e10
        )
        task.start()

    # Wait to not spam task submissions to GEE
    time.sleep(0.1)

def create_manifest(
        assetID,
        file_name_prefix,
        description_base,
        properties,
        bucket,
        year=None,
        model_start=None,
        model_end=None,
        i=None
):
    asset_name, image_name, _ = _create_name_description(
        assetID, file_name_prefix, description_base, year, i
    )

    if model_start is not None:
        start = model_start
        end = model_end + 1
    elif year is not None:
        start = year
        end = year + 1
    else:
        raise ValueError("Either year or model_start/model_end must be provided.")

    return {
        'name': asset_name,
        'properties': properties,
        'tilesets': [
            {'id': '0', 'sources': [ {'uris': [f'gs://{bucket}/{image_name}.tif']}]}
        ],
        'startTime': f'{start}-01-01T00:00:00.000000000Z',
        'endTime': f'{end}-01-01T00:00:00.000000000Z'
    }

def _create_name_description(assetID, file_name_prefix, description_base, year=None, i=None):
    suffix = ''
    if year is not None:
        suffix = f'_{year}'
    if i is not None:
        suffix += f'_tile_{i}'

    asset_name = assetID + suffix
    image_name = file_name_prefix + suffix
    description = description_base + suffix

    return asset_name, image_name, description