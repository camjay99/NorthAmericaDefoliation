import ee

def submit_job(
        image,
        assetID,
        file_name_prefix,
        description_base,
        year,
        scale,
        crs,
        region,
        cloudstorage,
        bucket=None,
        i=-1
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

def create_manifest(
        assetID,
        file_name_prefix,
        description_base,
        year,
        properties,
        bucket,
        i=-1
):
    asset_name, image_name, _ = _create_name_description(
        assetID, file_name_prefix, description_base, year, i
    )

    return {
        'name': asset_name,
        'properties': properties,
        'tilesets': [
            {'id': '0', 'sources': [ {'uris': [f'gs://{bucket}/{image_name}.tif']}]}
        ],
        'startTime': f'{year}-01-01T00:00:00.000000000Z',
        'endTime': f'{year+1}-01-01T00:00:00.000000000Z'
    }

def _create_name_description(assetID, file_name_prefix, description_base, year, i=-1):
    if i != -1:
        asset_name = f'{assetID}_{year}_tile_{i}'
        image_name = f'{file_name_prefix}_{year}_tile_{i}'
        description = f'{description_base}_{year}_tile_{i}'
    else:
        asset_name = f'{assetID}_{year}'
        image_name = f'{file_name_prefix}_{year}'
        description = f'{description_base}_{year}'

    return asset_name, image_name, description