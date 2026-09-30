import argparse
import os

import ee
import numpy as np
import rioxarray  # noqa: F401 - registers the .rio accessor
import xarray as xr
import xee  # noqa: F401 - registers the 'ee' xarray backend

import utils.geometries as geometries
import utils.preprocessing as preprocessing


##############################################################
# Parse arguments
##############################################################

parser = argparse.ArgumentParser(
    description='Options for calculating growing season phenology locally with xee.')

# The project to submit the code in.
# You may be prompted to to authenticate.
parser.add_argument('--project', '-p', action='store',
                    default=None, required=True)

# The first and last years to look for calculating average phenology.
parser.add_argument('--start', '-S', action='store', type=int, default=2019)
parser.add_argument('--end', '-E', action='store', type=int, default=2023)

# The data source to use for calculating trends.
parser.add_argument('--data', '-d', action='store',
                    default='HLS', choices=preprocessing.sources)

## NOTE: ONLY ONE OF GEOMETRY/STATE/RANGE CAN BE SPECIFIED
# The geometry to calculate trends within.
# A list of valid geometries are available in scripts/geometries.py
parser.add_argument('--geometry', '-g', action='store',
                    default=None, choices=geometries.site_names)

# State to calculate trends within.
parser.add_argument('--state', '-x', action='store', default=None)

# Use total Spongy Moth Range to calculate trends within.
parser.add_argument('--range', '-R', action='store_true')

# The crs to use for the output. Must be a projected crs with meter units
# since `scale` below (the HLS resolution) is specified in meters.
parser.add_argument('--crs', '-c', action='store', default='EPSG:5070')

# The width/length of grid cells to use for computation (in lat/lon degrees).
# Each grid cell is pulled from Earth Engine and processed locally, so this
# also controls how much data is held in memory at once.
parser.add_argument('--width', '-w', action='store', type=float, default=0.25)
parser.add_argument('--length', '-l', action='store', type=float, default=0.25)

# Directory to write local GeoTIFF output to.
parser.add_argument('--output_dir', '-o', action='store', default='.')

### Algorithm parameters
# Window size for calculating mean before/after conditions
parser.add_argument('--change_window', '-r', action='store', type=int, default=30)

# Threshold for initial binary classification
parser.add_argument('--threshold', '-t', action='store', type=float, default=0.5)

# Window size for smoothing data series
parser.add_argument('--smooth_window', '-W', action='store', type=int, default=3)

# Spacing between days of year evaluated as a potential change point.
# Increasing the value will make the algorithm run faster at the
# expense of lower fidelity in change dates.
parser.add_argument('--period', '-f', action='store', type=int, default=1)

# Parse arguments provided to script
args = parser.parse_args()


##############################################################
# Initialize Google Earth Engine API against the high-volume
# endpoint, which is tuned for many small/interactive requests
# like the per-pixel pulls xee makes below.
##############################################################

HIGH_VOLUME_ENDPOINT = 'https://earthengine-highvolume.googleapis.com'

try:
    ee.Initialize(project=args.project, opt_url=HIGH_VOLUME_ENDPOINT)
except Exception:
    # need to authenticate with your credential at the first time
    ee.Authenticate()
    ee.Initialize(project=args.project, opt_url=HIGH_VOLUME_ENDPOINT)


##################################################################
# Specify base names
##################################################################

num_not_specified = (args.geometry is None) + (args.state is None) + (not args.range)
assert num_not_specified != 3, "Must specify at least one of geometry, state, or range."
assert num_not_specified > 1, "Only one of geometry, state, and range can be specified at a time."

if args.geometry:
    name = args.geometry
    geometry = geometries.get_geometry(args.geometry)
elif args.state:
    name = args.state.replace(" ", "_")
    geometry = geometries.get_state(args.state)
else:
    name = "North_America"
    geometry = geometries.get_range()

scale = preprocessing.resolutions[args.data]
os.makedirs(args.output_dir, exist_ok=True)


##################################################################
# Split study regions into grid cells of specified size. Each
# cell is pulled from Earth Engine and processed locally in turn
# so that memory use stays bounded regardless of the size of
# `geometry`.
##################################################################

# Specify grid size in projection, x and y units (based on projection).
grid_projection = 'EPSG:4326'  # WGS84 lat lon
proj = ee.Projection(grid_projection).scale(args.width, args.length)
grid = geometry.coveringGrid(proj)

gridSize = grid.size().getInfo()
gridList = grid.toList(gridSize)

start_date = ee.Date.fromYMD(args.start, 1, 1)
end_date = ee.Date.fromYMD(args.end + 1, 1, 1)

doys = np.arange(0, 366, args.period)
change_steps = max(args.change_window, 1)


def circular_distance(a, b, period=365.0):
    """Shortest distance between two days-of-year on a 365-day circle."""
    diff = np.abs(a - b) % period
    return np.minimum(diff, period - diff)


def circular_forward_distance(a, b, period=365.0):
    """Distance travelled going forward in time from `a` to `b` (>= 0)."""
    return (b - a + period) % period


for i in range(gridSize):
    gridCell = ee.Feature(gridList.get(i)).geometry()

    ##################################################################
    # Fetch cloudmasked HLS imagery for this cell from Earth Engine
    # and pull it into a local, lazily-evaluated xarray Dataset via
    # xee/the high-volume endpoint. All cloud/water/shadow/forest-loss
    # masking happens server-side in preprocess_HLS; only the masked
    # EVI time series is transferred locally.
    ##################################################################

    if args.data == 'HLS':
        col = preprocessing.preprocess_HLS(
            start_date, end_date, gridCell, None, False, adddoy=False, aerosol_mask=3)

    # xee returns dims (time, y, x), with x/y coordinates in the requested
    # `crs`/`scale`.
    evi = xr.open_dataset(
        col.select('EVI'),
        engine='ee',
        crs=args.crs,
        scale=scale,
        geometry=gridCell,
    )

    if evi.sizes.get('time', 0) == 0:
        print(f"Tile {i}: no images found, skipping.")
        continue

    #evi = ds['EVI'].transpose('time', 'y', 'x')

    # (time) with values equal to day-of-year of the observation, in [0, 365).
    doy_of_obs = evi['time'].dt.dayofyear.values.astype(float) - 1.0

    ##################################################################
    # Window smoothing: for every candidate day-of-year, average all
    # (cloud-free) observations across all years that fall within
    # `smooth_window` days of it. This mirrors the per-doy
    # ee.Filter.dayOfYear + .mean() reduction in the GEE version, but
    # is computed once, locally, as a single weighted contraction
    # instead of one server round-trip per day.
    ##################################################################

    # Create array of shape (doy, doy_of_obs) using broadcasting. Entries in 
    # the matrix are 1 if the corresponding doy_of_obs is within `smooth_window` 
    # days of the corresponding doy, 0 otherwise. 
    smooth_weights = (circular_distance(doy_of_obs[None, :], doys[:, None])
                       <= args.smooth_window).astype(float)
    smooth_weights_da = xr.DataArray(
        smooth_weights, dims=['doy', 'time'],
        coords={'doy': doys, 'time': evi['time']})

    valid = evi.notnull().astype(float)
    evi_filled = evi.fillna(0.0)

    numerator = xr.dot(smooth_weights_da, evi_filled, dims='time')
    denom = xr.dot(smooth_weights_da, valid, dims='time')
    col_smoothed = (numerator / denom).where(denom > 0)

    ##################################################################
    # Estimate an absolute EVI threshold over the region from the
    # smoothed annual series, then convert to a binary "leaf-on" mask.
    ##################################################################

    max_EVI = col_smoothed.quantile(0.95, dim='doy', skipna=True)
    min_EVI = col_smoothed.quantile(0.05, dim='doy', skipna=True)
    amplitude = max_EVI - min_EVI
    thresh = amplitude * args.threshold + min_EVI

    col_bin = (col_smoothed > thresh).where(col_smoothed.notnull())

    ##################################################################
    # For every candidate day-of-year, compare the ratio of
    # "leaf-on" days in the `change_window` days before it to the
    # ratio in the `change_window` days after it.
    ##################################################################

    # M[d, d'] == 1 iff d' falls in (d, d + change_window], i.e. d' is
    # within `change_window` days after d (going forward, circularly).
    M = (circular_forward_distance(doys[:, None], doys[None, :])
         <= change_steps).astype(float)
    np.fill_diagonal(M, 0.0)

    after_weights = M
    # before_weights[d, d'] == 1 iff d' falls in [d - change_window, d)
    before_weights = M.T

    before_weights_da = xr.DataArray(
        before_weights, dims=['doy', 'doy_in'],
        coords={'doy': doys, 'doy_in': doys})
    after_weights_da = xr.DataArray(
        after_weights, dims=['doy', 'doy_in'],
        coords={'doy': doys, 'doy_in': doys})

    col_bin_valid = col_bin.notnull().astype(float)
    col_bin_filled = col_bin.fillna(0.0).astype(float)
    col_bin_in = col_bin_filled.rename({'doy': 'doy_in'})
    col_bin_valid_in = col_bin_valid.rename({'doy': 'doy_in'})

    # (n_doy, y, x) -> (n_doy, y, x), as it computes ratio for each doy.
    ratio_before = (xr.dot(before_weights_da, col_bin_in, dims='doy_in')
                     / xr.dot(before_weights_da, col_bin_valid_in, dims='doy_in'))
    ratio_after = (xr.dot(after_weights_da, col_bin_in, dims='doy_in')
                    / xr.dot(after_weights_da, col_bin_valid_in, dims='doy_in'))

    ratio_diff = ratio_before - ratio_after

    ##################################################################
    # SoS is the day-of-year where the ratio swings most sharply from
    # "leaf-off before" to "leaf-on after" (ratio_diff minimized).
    # EoS is the day-of-year where it swings the opposite way
    # (ratio_diff maximized).
    ##################################################################

    SoS = ratio_diff.idxmin(dim='doy', skipna=True).rename('SoS')
    EoS = ratio_diff.idxmax(dim='doy', skipna=True).rename('EoS')

    pheno = xr.Dataset({'SoS': SoS, 'EoS': EoS}).compute()
    pheno = pheno.fillna(0).astype('uint16')

    for var in ('SoS', 'EoS'):
        pheno[var].attrs = {}
    pheno.attrs.update({
        'method': 'Maximum Separability',
        'source': args.data,
        'start': args.start,
        'end': args.end,
        'project': 'NorthAmerica',
    })

    ##################################################################
    # Write results locally
    ##################################################################

    out_path = os.path.join(
        args.output_dir, f'{name}_Phenology_{args.data}_tile_{i}.tif')

    pheno_raster = pheno.to_array(dim='band')
    pheno_raster = pheno_raster.rio.set_spatial_dims(x_dim='x', y_dim='y')
    pheno_raster = pheno_raster.rio.write_crs(args.crs)
    pheno_raster.rio.to_raster(out_path)

    print(f"Tile {i}/{gridSize - 1}: wrote {out_path}")
