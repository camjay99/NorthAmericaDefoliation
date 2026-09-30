# The code in this module was generated with the assistance of
# Anthropic's Claude Sonnet 5 model. Code was then reviewed by
# Cameron Scholl to ensure it met the requirements of the project.

import math

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from rasterio.merge import merge
from rasterio.warp import transform_geom, transform_bounds
from rasterio.enums import Resampling
from rasterio.features import geometry_mask
from rasterio.windows import Window, bounds as window_bounds
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


def compute_grid(polygon, polygon_crs, dst_crs, resolution):
    reprojected = shape(transform_geom(polygon_crs, dst_crs, mapping(polygon)))
    minx, miny, maxx, maxy = reprojected.bounds

    width = max(1, math.ceil((maxx - minx) / resolution))
    height = max(1, math.ceil((maxy - miny) / resolution))
    transform = from_origin(minx, maxy, resolution, resolution)

    return transform, width, height


def get_raster_index(filepaths, dst_crs):
    assert len(filepaths) > 0, "filepaths must not be empty"

    index = []
    for filepath in filepaths:
        with rasterio.open(filepath) as dataset:
            bounds = transform_bounds(dataset.crs, dst_crs, *dataset.bounds)
        index.append((filepath, bounds))

    return index


def iter_windows(width, height, block_size):
    for row_off in range(0, height, block_size):
        block_height = min(block_size, height - row_off)
        for col_off in range(0, width, block_size):
            block_width = min(block_size, width - col_off)
            yield Window(col_off, row_off, block_width, block_height)


def sample_polygons(filepath, polygons, polygons_crs, bands=None,
                     mask_path=None, mask_band=1, mask_valid_value=1):
    # Returns a DataFrame with columns ['polygon', 'row', 'col',
    # 'band_{x}', ...] where each row corresponds to a pixel in the raster that
    # is fully contained within one of the polygons.
    #
    # If mask_path is given (e.g. a Planet UDM2 file), it is assumed to
    # share the same pixel grid as filepath (true for a UDM2 delivered
    # alongside its analytic image). Pixels whose mask_band value is not
    # equal to mask_valid_value (default: band 1 "clear" == 1) are
    # excluded, e.g. to avoid sampling invariant targets under cloud/
    # shadow/haze/snow.
    assert len(polygons) > 0, "polygons must not be empty"

    with rasterio.open(filepath) as dataset:
        band_indexes = bands if bands is not None else list(range(1, dataset.count + 1))
        nodata = dataset.nodata
        raster_bounds = box(*dataset.bounds)

        records = []
        for polygon_idx, polygon in enumerate(polygons):
            reprojected = shape(transform_geom(polygons_crs, dataset.crs, mapping(polygon)))
            if reprojected.is_empty or not reprojected.intersects(raster_bounds):
                continue

            # Read a window covering the polygon (padded by a pixel on each
            # side), then keep only pixels whose full cell - not just their
            # center - falls within the polygon.
            minx, miny, maxx, maxy = reprojected.bounds
            row_start, col_start = dataset.index(minx, maxy, op=math.floor)
            row_stop, col_stop = dataset.index(maxx, miny, op=math.ceil)
            row_start, col_start = row_start - 1, col_start - 1
            row_stop, col_stop = row_stop + 1, col_stop + 1

            window = Window(col_start, row_start,
                            col_stop - col_start, row_stop - row_start)
            window_transform = dataset.window_transform(window)
            out_shape = (int(window.height), int(window.width))

            # A pixel only counts as "in" the polygon if its whole cell is
            # inside it. Build the inverse of the polygon (the window's
            # extent with the polygon cut out of it) and rasterize that
            # with all_touched=True: any pixel touching that inverse
            # region is at least partially outside the polygon, so the
            # remaining (untouched) pixels are exactly the ones fully
            # contained by the polygon.
            inverse_polygon = box(*window_bounds(window, dataset.transform)).difference(reprojected)
            if inverse_polygon.is_empty:
                touches_outside = np.zeros(out_shape, dtype=bool)
            else:
                touches_outside = geometry_mask(
                    [inverse_polygon], out_shape=out_shape,
                    transform=window_transform, all_touched=True, invert=True)
            contained = ~touches_outside
            if not contained.any():
                continue

            local_rows, local_cols = np.nonzero(contained)

            block = dataset.read(band_indexes, window=window,
                                    boundless=True, out_dtype='float64')
            values = block[band_indexes, local_rows, local_cols]
            valid = np.all(np.isfinite(values), axis=0)
            if nodata is not None:
                valid &= values != nodata

            if mask_dataset is not None:
                with rasterio.open(mask_path) as mask_dataset:
                    mask_values = mask_dataset.read(
                        mask_band, window=window, boundless=True, fill_value=0)
                    valid &= mask_values[local_rows, local_cols] == mask_valid_value

            if not valid.any():
                continue

            info = {
                'polygon': polygon_idx,
                'row': local_rows[valid],
                'col': local_cols[valid],
                }
            bands_out = {f'band_{band_idx}': values[valid, i]
                            for i, band_idx in enumerate(band_indexes)}
            records.append(pd.DataFrame(info | bands_out))

        if not records:
            bands_columns = [f'band_{band_idx}' for band_idx in band_indexes]
            return pd.DataFrame(columns=['polygon', 'row', 'col'] + bands_columns)

        return pd.concat(records, ignore_index=True)