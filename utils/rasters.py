# The code in this module was generated with the assistance of 
# Anthropic's Claude Sonnet 5 model. Code was then reviewed by 
# Cameron Scholl to ensure it met the requirements of the project.

import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.merge import merge
from rasterio.warp import transform_geom
from rasterio.enums import Resampling
from shapely.geometry import box, mapping, shape


def filter_intersecting_rasters(filepaths, polygon, polygon_crs):
    assert len(filepaths) > 0, "filepaths must not be empty"

    intersecting_filepaths = []
    for filepath in filepaths:
        with rasterio.open(filepath) as dataset:
            reprojected_polygon = shape(transform_geom(polygon_crs, dataset.crs, mapping(polygon)))
            if reprojected_polygon.intersects(box(*dataset.bounds)):
                intersecting_filepaths.append(filepath)

    return intersecting_filepaths


def open_rasters_as_vrt(filepaths, dst_crs, **kwargs):
    assert len(filepaths) > 0, "filepaths must not be empty"

    vrts = []
    for filepath in filepaths:
        dataset = rasterio.open(filepath)
        vrts.append(WarpedVRT(dataset, crs=dst_crs, **kwargs))

    return vrts


def mosaic_rasters(datasets, polygon, polygon_crs, method='first',
                   height=None, width=None, resampling=Resampling.nearest):
    assert len(datasets) > 0, "datasets must not be empty"
    assert (height is None) == (width is None), "height and width must be given together"

    reprojected_polygon = shape(transform_geom(polygon_crs, datasets[0].crs, mapping(polygon)))
    bounds = reprojected_polygon.bounds

    res = None
    if height is not None:
        res = ((bounds[2] - bounds[0]) / width, (bounds[3] - bounds[1]) / height)

    mosaic_array, mosaic_transform = merge(datasets, bounds=bounds, res=res,
                                            resampling=resampling, method=method)

    return mosaic_array, mosaic_transform
