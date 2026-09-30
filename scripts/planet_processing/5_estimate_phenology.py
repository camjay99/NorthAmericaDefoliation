# The code in this module was generated with the assistance of
# Anthropic's Claude Sonnet 5 model based on the maximum separation
# algorithm for Google Earth Engine presented in Descals et al. 2020. 
# Code was then reviewed by Cameron Scholl to ensure it met the 
# requirements of the project.

import argparse
import datetime
import glob
import os
import re

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import bounds as window_bounds

import utils.rasters as rasters

NODATA = 65535

# Matches the {year}{month}{day}_ prefix used for PlanetScope scene filenames.
DATE_PATTERN = re.compile(r'^(\d{4})(\d{2})(\d{2})_')


##############################################################
# Parse arguments
##############################################################

parser = argparse.ArgumentParser(
    description='Estimate average start/end of season (SoS/EoS) locally from '
                'PlanetScope imagery within a specified polygon using '
                'the maximum separability method.')

parser.add_argument('--input', '-i', action='store', required=True,
                     help='Directory containing PlanetScope analytic '
                          'imagery. Filenames must start with a date in '
                          '{year}{month}{day} form, e.g. 20230615_*.tif')

parser.add_argument('--pattern', action='store', default='*.tif',
                     help='Glob pattern (relative to --input) selecting '
                          'which files to use, e.g. "*_AnalyticMS.tif" if '
                          'the folder also has UDM or metadata files.')

parser.add_argument('--geometry', '-g', action='store', required=True,
                     help='Path to a vector file (e.g. GeoJSON/shapefile) '
                          'giving the polygon to estimate phenology within.')

parser.add_argument('--output', '-o', action='store', required=True,
                     help='Path to save the output SoS/EoS GeoTIFF.')

# The first and last years of imagery to include when averaging phenology.
parser.add_argument('--start', '-S', action='store', type=int, default=None)
parser.add_argument('--end', '-E', action='store', type=int, default=None)

# The crs/resolution to use for the output.
parser.add_argument('--crs', '-c', action='store', default='EPSG:5070')
parser.add_argument('--resolution', '-res', action='store', type=float,
                     default=3.0, help='Output pixel size, in --crs units.')

# Band indexes (1-based) for a 4-band PlanetScope analytic product.
parser.add_argument('--blue-band', action='store', type=int, default=1)
parser.add_argument('--red-band', action='store', type=int, default=3)
parser.add_argument('--nir-band', action='store', type=int, default=4)

# UDM2 (usable data mask) cloud/shadow/haze/snow masking. Only pixels
# flagged "clear" (band 1 == 1) are used in the calculations.
parser.add_argument('--udm-dir', action='store', default=None,
                     help='Directory containing UDM2 files, if different '
                          'from --input.')
parser.add_argument('--udm-search', action='store', default=r'AnalyticMS(_SR)?\.tif$',
                     help='Regex matched against each analytic filename '
                          'and replaced with --udm-replace to find its '
                          'matching UDM2 file.')
parser.add_argument('--udm-replace', action='store', default='udm2.tif',
                     help='Replacement text used with --udm-search to '
                          'build the UDM2 filename.')

# The width/height (in pixels) of blocks used for computation, so that all
# imagery needed to fit a block's model fits comfortably in memory.
parser.add_argument('--block-size', action='store', type=int, default=256)

### Algorithm parameters
# Window size for calculating mean before/after conditions
parser.add_argument('--change_window', '-r', action='store', type=int, default=30)

# Threshold for initial binary classification
parser.add_argument('--threshold', '-t', action='store', type=float, default=0.5)

# Window size for smoothing data series
parser.add_argument('--smooth_window', '-W', action='store', type=int, default=3)

# Spacing between days of year evaluated as a potential change point.
# Increasing the value will make the algorithm run faster at the
# expense of lower fidelity in change dates
parser.add_argument('--period', '-f', action='store', type=int, default=1)

args = parser.parse_args()

assert 0 < args.threshold < 1, "--threshold must be between 0 and 1."


##############################################################
# Load region of interest
##############################################################

roi = gpd.read_file(args.geometry)
polygon = roi.union_all()
polygon_crs = roi.crs


##############################################################
# Find input files and extract day-of-year metadata
##############################################################

image_filepaths = sorted(glob.glob(os.path.join(args.input, args.pattern)))
assert image_filepaths, f'No files matching {args.pattern} found in {args.input}.'

image_filepaths = rasters.filter_intersecting_rasters(image_filepaths, polygon, polygon_crs)
assert image_filepaths, 'No input files intersect --geometry.'

doy_by_path = {}
for filepath in image_filepaths:
    match = DATE_PATTERN.match(os.path.basename(filepath))
    assert match, (f'{filepath} does not start with a {{year}}{{month}}'
                    '{day}_ date prefix.')
    year, month, day = (int(g) for g in match.groups())
    if args.start is not None and year < args.start:
        continue
    if args.end is not None and year > args.end:
        continue
    date = datetime.date(year, month, day)
    doy_by_path[filepath] = (date - datetime.date(year, 1, 1)).days

assert doy_by_path, 'No input files remain after applying --start/--end.'
print(f'Using {len(doy_by_path)} image(s) spanning day-of-year '
      f'{min(doy_by_path.values())}-{max(doy_by_path.values())}.')


##############################################################
# Locate each image's matching UDM2 cloud/shadow/haze/snow mask
##############################################################

udm_dir = args.udm_dir or args.input
udm_by_path = {}
for filepath in doy_by_path:
    basename = os.path.basename(filepath)
    udm_name, n_subs = re.subn(args.udm_search, args.udm_replace, basename)
    assert n_subs > 0, (
        f'--udm-search ({args.udm_search!r}) did not match {basename!r}. '
        'Pass --udm-search/--udm-replace matching your file naming.')
    udm_path = os.path.join(udm_dir, udm_name)
    assert os.path.isfile(udm_path), (
        f'No UDM2 file found for {filepath} (expected {udm_path}).')
    udm_by_path[filepath] = udm_path


##############################################################
# Build output grid and pre-index raster bounds
##############################################################

transform, width, height = rasters.compute_grid(
    polygon, polygon_crs, args.crs, args.resolution)

bounds_index = rasters.get_raster_index(list(doy_by_path), args.crs)


##############################################################
# Precompute the day-of-year grid and change-window matrices
##############################################################

def circular_window(values, start, end):
    # Boolean mask for `values` (day-of-year) within [start, end] mod 365,
    # wrapping around the end of the year when start > end. 
    values = np.asarray(values)
    start = start % 365
    end = end % 365
    if start <= end:
        return (values >= start) & (values <= end)
    else:
        return (values >= start) | (values <= end)


list_dates = np.arange(0, 366, args.period, dtype=np.float64)
n_dates = len(list_dates)

# (n_dates, n_dates) membership matrices for the before/after windows used
# to estimate the ratio of "green" observations around each candidate day.
# These are determined by the --change_window parameter, and will not change
# accross blocks.
before_matrix = np.zeros((n_dates, n_dates), dtype=np.float32)
after_matrix = np.zeros((n_dates, n_dates), dtype=np.float32)
for i, day in enumerate(list_dates):
    before_matrix[i] = circular_window(list_dates, day - args.change_window, day)
    after_matrix[i] = circular_window(list_dates, day, day + args.change_window)


##############################################################
# Per-block computation
##############################################################

def bounds_intersect(a, b):
    return not (a[2] <= b[0] or a[0] >= b[2] or a[3] <= b[1] or a[1] >= b[3])


def process_block(block_files, window):
    block_shape = (int(window.height), int(window.width))
    if not block_files:
        empty = np.full(block_shape, NODATA, dtype=np.uint16)
        return empty, empty.copy()

    evi_stack = []
    valid_stack = []
    doys = []
    for filepath, doy in block_files:
        with rasterio.open(filepath) as src:
            with WarpedVRT(src, crs=args.crs, transform=transform,
                            width=width, height=height,
                            resampling=Resampling.bilinear) as vrt:
                blue, red, nir = vrt.read(
                    [args.blue_band, args.red_band, args.nir_band],
                    window=window, out_dtype='float32',
                    boundless=True, fill_value=np.nan)
                nodata = vrt.nodata

        valid = np.isfinite(blue) & np.isfinite(red) & np.isfinite(nir)
        if nodata is not None:
            valid &= (blue != nodata) & (red != nodata) & (nir != nodata)

        udm_path = udm_by_path[filepath]
        with rasterio.open(udm_path) as src:
            with WarpedVRT(src, crs=args.crs, transform=transform,
                            width=width, height=height,
                            resampling=Resampling.nearest) as vrt:
                clear = vrt.read(1, window=window, boundless=True,
                                  fill_value=0)
        valid &= (clear == 1)

        with np.errstate(divide='ignore', invalid='ignore'):
            evi = 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
        valid &= np.isfinite(evi) & (evi >= 0) & (evi <= 1)

        evi_stack.append(np.where(valid, evi, 0.0))
        valid_stack.append(valid)
        doys.append(doy)

    # n_images x (height, width) -> (n_images, height, width)
    evi_stack = np.stack(evi_stack, axis=0)
    valid_stack = np.stack(valid_stack, axis=0).astype(np.float32)
    doys = np.array(doys, dtype=np.float64)

    #################################
    # Window smoothing
    #################################
    smooth_matrix = np.zeros((n_dates, len(doys)), dtype=np.float32)
    for i, day in enumerate(list_dates):
        smooth_matrix[i] = circular_window(
            doys, day - args.smooth_window, day + args.smooth_window)

    # d = n_dates, n = n_images, h = height, w = width
    # We sum accross the n_images axis in both arrays. Smooth_matrix is 1
    # if a given image is in the smoothing window for a given day, 
    # and 0 otherwise. Thus, this is a moving window sum.
    smoothed_sum = np.einsum('dn,nhw->dhw', smooth_matrix, evi_stack)
    smoothed_count = np.einsum('dn,nhw->dhw', smooth_matrix, valid_stack)
    with np.errstate(divide='ignore', invalid='ignore'):
        smoothed = smoothed_sum / smoothed_count
    smoothed_valid = smoothed_count > 0

    ############################################
    # Estimate absolute threshold over the block
    ############################################
    smoothed_nan = np.where(smoothed_valid, smoothed, np.nan)
    max_evi = np.nanpercentile(smoothed_nan, 95, axis=0)
    min_evi = np.nanpercentile(smoothed_nan, 5, axis=0)
    thresh = (max_evi - min_evi) * args.threshold + min_evi

    valid_f = smoothed_valid.astype(np.float32)
    binary = np.where(smoothed_valid & (smoothed > thresh[np.newaxis]), 1.0, 0.0)
    binary = binary.astype(np.float32)

    ##################################################
    # Estimate ratio of observations above the
    # threshold before and after each day of the year
    ##################################################
    # d = n_dates, e = n_dates, h = height, w = width
    # In before_matrix/after_matrix, entry (i, j) is 1 if day j is in 
    # the before/after window for day i, and 0 otherwise.
    before_sum = np.einsum('de,ehw->dhw', before_matrix, binary)
    before_count = np.einsum('de,ehw->dhw', before_matrix, valid_f)
    after_sum = np.einsum('de,ehw->dhw', after_matrix, binary)
    after_count = np.einsum('de,ehw->dhw', after_matrix, valid_f)

    with np.errstate(divide='ignore', invalid='ignore'):
        mean_before = before_sum / before_count
        mean_after = after_sum / after_count
    mean_diff = mean_before - mean_after

    ##################################################
    # Extract SoS and EoS
    ##################################################
    # Ensure that no data dates can never be picked.
    has_data = ~np.all(np.isnan(mean_diff), axis=0)
    diff_min = np.where(np.isnan(mean_diff), np.inf, mean_diff)
    diff_max = np.where(np.isnan(mean_diff), -np.inf, mean_diff)
    sos_idx = np.argmin(diff_min, axis=0)
    eos_idx = np.argmax(diff_max, axis=0)

    sos = np.where(has_data, list_dates[sos_idx], NODATA).astype(np.uint16)
    eos = np.where(has_data, list_dates[eos_idx], NODATA).astype(np.uint16)

    return sos, eos


##############################################################
# Process the output raster in blocks and write results
##############################################################

profile = {
    'driver': 'GTiff',
    'dtype': 'uint16',
    'nodata': NODATA,
    'width': width,
    'height': height,
    'count': 2,
    'crs': args.crs,
    'transform': transform,
    'tiled': True,
    'blockxsize': 256,
    'blockysize': 256,
    'compress': 'deflate',
}

output_dir = os.path.dirname(os.path.abspath(args.output))
os.makedirs(output_dir, exist_ok=True)

with rasterio.open(args.output, 'w', **profile) as dst:
    dst.set_band_description(1, 'SoS')
    dst.set_band_description(2, 'EoS')
    dst.update_tags(method='Maximum Separability', source='PlanetScope',
                     start=args.start, end=args.end)

    windows = list(rasters.iter_windows(width, height, args.block_size))
    for block_num, window in enumerate(windows, 1):
        block_extent = window_bounds(window, transform)
        block_files = [
            (filepath, doy_by_path[filepath])
            for filepath, bounds in bounds_index
            if bounds_intersect(bounds, block_extent)
        ]

        sos_block, eos_block = process_block(block_files, window)

        dst.write(sos_block, 1, window=window)
        dst.write(eos_block, 2, window=window)
        print(f'Processed block {block_num}/{len(windows)} '
              f'({len(block_files)} image(s) contributing).')

print(f'Saved phenology estimate to {args.output}.')
