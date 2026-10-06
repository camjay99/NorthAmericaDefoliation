# The code in this module was generated with the assistance of
# Anthropic's Claude Sonnet 5 model. Code was then reviewed by
# Cameron Scholl to ensure it met the requirements of the project.

import argparse
import glob
import os
import re

import geopandas as gpd
import numpy as np
import rasterio

import utils.rasters as rasters


##############################################################
# Parse arguments
##############################################################

parser = argparse.ArgumentParser(
    description='Filter out poorly georegistered PlanetScope images and '
                 'apply pseudo-invariant-feature radiometric correction to '
                 'the remainder, relative to a high quality reference '
                 'image.')

parser.add_argument('--input', '-i', action='store', required=True,
                     help='Directory containing candidate images to '
                          'filter/correct.')

parser.add_argument('--pattern', action='store', default='*.tif',
                     help='Glob pattern (relative to --input) selecting '
                          'which files to consider.')

parser.add_argument('--reference', '-r', action='store', required=True,
                     help='Path to a high quality, well georegistered and '
                          'radiometrically trusted image used as the basis '
                          'for evaluating/correcting the other images.')

parser.add_argument('--geo-targets', action='store', required=True,
                     help='Path to a vector file of invariant target '
                          'polygons (e.g. bright, stable features) used to '
                          'evaluate georegistration.')

parser.add_argument('--radio-targets', action='store', required=True,
                     help='Path to a vector file of pseudo-invariant target '
                          'polygons used to fit the radiometric '
                          'correction.')

parser.add_argument('--output', '-o', action='store', required=True,
                     help='Directory to save the filtered/corrected images '
                          'to.')

parser.add_argument('--geo-r2-threshold', action='store', type=float,
                     default=0.8, help='Minimum R2, fit between an image '
                          "and --reference's NIR values at --geo-targets, "
                          'required to keep an image.')

parser.add_argument('--radio-r2-threshold', action='store', type=float,
                     default=0.7, help='Minimum R2, fit per-band between an '
                          "image and --reference's values at "
                          '--radio-targets, required to keep an image.')

parser.add_argument('--nir-band', action='store', type=int, default=4,
                     help='1-based band index of the NIR band, used to '
                          'evaluate georegistration.')

parser.add_argument('--min-points', action='store', type=int, default=3,
                     help='Minimum number of valid (fully-contained, '
                          'non-nodata) image pixels required to fit a '
                          'regression.')

parser.add_argument('--block-size', action='store', type=int, default=1024,
                     help='Width/height, in pixels, of the blocks used '
                          'when writing corrected imagery.')

# UDM2 (usable data mask) cloud/shadow/haze/snow masking. Only pixels
# flagged "clear" (band 1 == 1) are used when sampling --geo-targets and
# --radio-targets, in case clouds happen to cover an invariant target.
parser.add_argument('--udm-dir', action='store', default=None,
                     help='Directory containing UDM2 files matching '
                          '--input images and --reference, if different '
                          'from --input. If omitted, no UDM masking is '
                          'applied.')
parser.add_argument('--udm-search', action='store', default=r'AnalyticMS(_SR)?\.tif$',
                     help='Regex matched against each analytic filename '
                          'and replaced with --udm-replace to find its '
                          'matching UDM2 file.')
parser.add_argument('--udm-replace', action='store', default='udm2.tif',
                     help='Replacement text used with --udm-search to '
                          'build the UDM2 filename.')

args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)


##############################################################
# Load invariant target polygons
##############################################################

geo_targets = gpd.read_file(args.geo_targets)
assert len(geo_targets) > 0, '--geo-targets must contain at least one polygon.'

radio_targets = gpd.read_file(args.radio_targets)
assert len(radio_targets) > 0, '--radio-targets must contain at least one polygon.'


##############################################################
# Find candidate images
##############################################################

image_filepaths = sorted(glob.glob(os.path.join(args.input, args.pattern)))
reference_path = os.path.abspath(args.reference)
image_filepaths = [
    filepath for filepath in image_filepaths
    if os.path.abspath(filepath) != reference_path
]
assert image_filepaths, f'No files matching {args.pattern} found in {args.input}.'
print(f'Found {len(image_filepaths)} candidate image(s).')

with rasterio.open(args.reference) as dataset:
    reference_band_count = dataset.count


##############################################################
# Locate each image's matching UDM2 cloud/shadow/haze/snow mask
##############################################################

udm_by_path = {}
if args.udm_dir is not None:
    udm_dir = args.udm_dir
    for filepath in image_filepaths:
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
# Fit a linear regression between two sets of sampled values
##############################################################

def fit_line(x, y, min_points):
    # Fits both an ordinary least-squares R2 (used to evaluate goodness of
    # fit) and a reduced major axis (RMA) slope/intercept (used to correct
    # x onto y's scale, appropriate since both x and y are noisy
    # measurements rather than one being an error-free predictor).
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    if x.size < min_points:
        return None

    std_x = x.std(ddof=1)
    std_y = y.std(ddof=1)
    if std_x == 0 or std_y == 0:
        return None

    r = np.corrcoef(x, y)[0, 1]
    if not np.isfinite(r):
        return None
    r2 = r ** 2

    slope = np.sign(r) * (std_y / std_x)
    intercept = y.mean() - slope * x.mean()

    return slope, intercept, r2


##############################################################
# Evaluate georegistration and fit radiometric correction
##############################################################

JOIN_KEYS = ['polygon', 'row', 'col']

geo_polygons = geo_targets.geometry.tolist()
reference_geo_df = rasters.sample_polygons(
    args.reference, geo_polygons, geo_targets.crs, bands=[args.nir_band])

radio_polygons = radio_targets.geometry.tolist()
reference_radio_df = rasters.sample_polygons(
    args.reference, radio_polygons, radio_targets.crs)

kept_filepaths = []
corrections_by_path = {}

for filepath in image_filepaths:
    with rasterio.open(filepath) as dataset:
        band_count = dataset.count
    if band_count != reference_band_count:
        print(f'Removed {filepath}: has {band_count} band(s), expected '
              f'{reference_band_count} to match --reference.')
        continue

    #################################
    # Evaluate georegistration
    #################################
    image_geo_df = rasters.sample_polygons(
        filepath, geo_polygons, geo_targets.crs, bands=[args.nir_band],
        mask_path=udm_by_path.get(filepath))
    geo_pairs = image_geo_df.merge(reference_geo_df, on=JOIN_KEYS,
                                    suffixes=('_image', '_reference'))

    geo_fit = fit_line(geo_pairs[f'band_{args.nir_band}_image'].to_numpy().astype(float),
                        geo_pairs[f'band_{args.nir_band}_reference'].to_numpy().astype(float),
                        args.min_points)
    if geo_fit is None:
        print(f'Removed {filepath}: too few valid --geo-targets pixels.')
        continue

    _, _, geo_r2 = geo_fit
    if geo_r2 < args.geo_r2_threshold:
        print(f'Removed {filepath}: georegistration R2 {geo_r2:.3f} < '
              f'{args.geo_r2_threshold}.')
        continue

    #################################
    # Fit per-band radiometric correction
    #################################
    image_radio_df = rasters.sample_polygons(
        filepath, radio_polygons, radio_targets.crs,
        mask_path=udm_by_path.get(filepath))
    radio_pairs = image_radio_df.merge(reference_radio_df, on=JOIN_KEYS,
                                        suffixes=('_image', '_reference'))

    band_corrections = []
    band_r2s = []
    failed_band = None
    for band_idx in range(1, band_count + 1):
        radio_fit = fit_line(radio_pairs[f'band_{band_idx}_image'].to_numpy().astype(float),
                              radio_pairs[f'band_{band_idx}_reference'].to_numpy().astype(float),
                              args.min_points)
        if radio_fit is None:
            failed_band = band_idx
            break

        slope, intercept, radio_r2 = radio_fit
        if radio_r2 < args.radio_r2_threshold:
            failed_band = band_idx
            break

        band_corrections.append((slope, intercept))
        band_r2s.append(radio_r2)

    if failed_band is not None:
        print(f'Removed {filepath}: radiometric correction failed for '
              f'band {failed_band} (R2 below {args.radio_r2_threshold} '
              'or too few valid points).')
        continue

    band_r2_summary = ', '.join(
        f'band {i + 1} R2 {r2:.3f}' for i, r2 in enumerate(band_r2s))
    print(f'Kept {filepath}: georegistration R2 {geo_r2:.3f}, '
          f'{band_r2_summary}.')
    kept_filepaths.append(filepath)
    corrections_by_path[filepath] = band_corrections


##############################################################
# Apply the radiometric correction and save the corrected images
##############################################################

for filepath in kept_filepaths:
    band_corrections = corrections_by_path[filepath]
    out_path = os.path.join(args.output, os.path.basename(filepath))

    with rasterio.open(filepath) as src:
        profile = src.profile.copy()
        nodata = src.nodata

        # Range to clip corrected values to
        dtype = np.dtype(profile['dtype'])
        if np.issubdtype(dtype, np.integer):
            info = np.iinfo(dtype)
            value_min, value_max = info.min, info.max
        else:
            value_min, value_max = -np.inf, np.inf

        with rasterio.open(out_path, 'w', **profile) as dst:
            windows = list(rasters.iter_windows(
                src.width, src.height, args.block_size))
            for window in windows:
                block = src.read(window=window)
                valid = np.ones(block.shape[1:], dtype=bool)
                if nodata is not None:
                    valid &= np.all(block != nodata, axis=0)

                corrected = block.astype('float64')
                for band_idx, (slope, intercept) in enumerate(band_corrections):
                    corrected[band_idx] = corrected[band_idx] * slope + intercept

                corrected = np.clip(corrected, value_min, value_max)
                corrected = corrected.astype(dtype)
                if nodata is not None:
                    corrected = np.where(valid[np.newaxis], corrected, nodata)

                dst.write(corrected, window=window)

    print(f'Saved corrected image to {out_path}.')

print(f'Kept {len(kept_filepaths)}/{len(image_filepaths)} image(s).')