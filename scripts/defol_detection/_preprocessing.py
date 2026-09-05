import ee

sources = ['HLS']
resolutions = {'HLS':30}    


def preprocess_HLS(start_date, end_date, geometry, 
                   start_doy=0, end_doy=365, phenology=None, adddoy=True):
    # Load HLS S30/L30
    collection_L30 = (ee.ImageCollection("NASA/HLS/HLSL30/v002")
                        .filterDate(start_date, end_date)
                        .filter(ee.Filter.dayOfYear(start_doy, end_doy))
                        .filterBounds(geometry))
                    
    collection_S30 = (ee.ImageCollection("NASA/HLS/HLSS30/v002")
                        .filterDate(start_date, end_date)
                        .filter(ee.Filter.dayOfYear(start_doy, end_doy))
                        .filterBounds(geometry))
    
    # Link the Cloud Score + collection to the S30 collection
    cloud_score = (ee.ImageCollection("GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED")
                    .filterDate(start_date, end_date)
                    .filter(ee.Filter.dayOfYear(start_doy, end_doy))
                    .filterBounds(geometry))
    
    def reindex(image):
        old_index = ee.String(image.get('system:index'))
        tags = old_index.split('_')
        tile_tag = ee.String(tags.get(2))
        date_tag = ee.String(tags.get(0))
        new_index = tile_tag.cat('_').cat(date_tag)
        return image.set('HLS_index', new_index)
    cloud_score = cloud_score.map(reindex)

    filter = ee.Filter.equals(
        leftField='system:index',
        rightField='HLS_index')    
    join = ee.Join.inner()
    joined = join.apply(collection_S30, cloud_score, filter)
    collection_S30 = ee.ImageCollection(
        joined.map(
            lambda feature: ee.Image.cat(feature.get('primary'), 
                                        feature.get('secondary'))))

    # Rename bands so collections can be used together.
    def rename_bands_L30(image):
        return image.select(['B2','B3','B4','B5', 'B6', 'B7', 'Fmask'],
                            ['BLUE','GREEN','RED','NIR', 'SWIR1', 'SWIR2', 'Fmask'])

    def rename_bands_S30(image):
        return image.select(['B2','B3','B4','B8', 'B11', 'B12', 'Fmask','cs'],
                            ['BLUE','GREEN','RED','NIR','SWIR1','SWIR2','Fmask','cs'])
    
    # Load Hansen et al. 2013 dataset for masking forest change
    forest_change = ee.Image("UMD/hansen/global_forest_change_2025_v1_13")
        
    def preprocess(image, s2=False):
        # New Bands
        EVI = image.expression(
            '2.5 * ((NIR-RED) / (NIR + 6 * RED - 7.5* BLUE + 1))', {
                'NIR': image.select('NIR'),
                'RED': image.select('RED'),
                'BLUE': image.select('BLUE')
            }).rename('EVI')
        doy = image.date().getRelative('day', 'year')
        if adddoy or (phenology != None):
            doy_band = ee.Image.constant(doy).uint16().rename('doy')
        year = image.date().get('year').subtract(2000) # Global Forest Change uses 00-24
        
        ### Masks
        ## Fmask cloud mask
        # Bit 0 - Cirrus (unused)
        # Bit 1 - Cloud
        # Bit 2 - Adjacent to cloud/shadow
        # Bit 3 - Cloud shadow
        # Bit 4 - Snow/Ice
        # Bit 5 - Water
        # Bit 6-7 Aerosol level
        # We mask clouds, shadows, adjecent areas, water and moderate/high aerosol.
        # Note that small clouds/shadows are often incorrectly classified as aerosols.
        fmask_val = ee.Number.parse("00101110", 2)
        mask = image.select("Fmask").bitwiseAnd(fmask_val).eq(0) 
        ## Cloud Score+ for HLS S30 images
        # Fmask fails to detect cloud shadows and haze quite frequently for 
        # Sentinel-2 images, thus we also rely on a dedicated usability score 
        # that captures most of these failure points.
        if s2:
            mask = mask.And(image.select('cs').gte(0.60))
        #albedo = image.select('BLUE').add(image.select('RED')).add(image.select('GREEN')) Maybe unneeded with Cloud Score+
        #mask = mask.And(albedo.lte(0.75)) # Albedo mask in case Fmask failed
        mask = mask.And(EVI.lte(1).And(EVI.gte(0))) # EVI mask
        mask = mask.And(
            forest_change.select('lossyear')
                         .lte(year).unmask().Not()) # Forest change mask
        if phenology != None:
            mask = mask.And(doy_band.gte(phenology.select('SoS')) # Phenology mask
                            .And(doy_band.lte(phenology.select('EoS'))))

        SWIR = image.select(['SWIR1', 'SWIR2'])
        if adddoy:
            result = ee.Image([EVI, SWIR, doy_band])
        else:
            result = ee.Image([EVI, SWIR])
        return (result.updateMask(mask)
                     .copyProperties(image, ['system:time_start']))

    collection_L30 = collection_L30.map(
        lambda image: preprocess(rename_bands_L30(image)))
    collection_S30 = collection_S30.map(
        lambda image: preprocess(rename_bands_S30(image), True))

    return collection_L30.merge(collection_S30)