# The code in this module was generated with the assistance of 
# Anthropic's Claude Sonnet 5 model. Code was then reviewed by 
# Cameron Scholl to ensure it met the requirements of the project.

import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.merge import merge
from rasterio.windows import from_bounds
from rasterio.windows import transform as window_transform
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


def mosaic_rasters(datasets, polygon,
                   polygon_crs, method='first'):
    assert len(datasets) > 0, "datasets must not be empty"

    reprojected_polygon = shape(transform_geom(polygon_crs, datasets[0].crs, mapping(polygon)))
    mosaic_array, mosaic_transform = merge(datasets, method=method)

    window = from_bounds(*reprojected_polygon.bounds, transform=mosaic_transform).round_lengths().round_offsets()
    row_slice, col_slice = window.toslices()
    clipped_array = mosaic_array[:, row_slice, col_slice]
    clipped_transform = window_transform(window, mosaic_transform)

    return clipped_array, clipped_transform
