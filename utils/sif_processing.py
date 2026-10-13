import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio as rio
import shapely
from scipy.stats import theilslopes

def sample_raster(gdf, raster_path):
    """
    Sample a raster file using geometries from a GeoDataFrame and return a 
    filtered GeoDataFrame containing only the geometries that intersect with the 
    raster.

    Parameters:
    gdf (gpd.GeoDataFrame): Input GeoDataFrame to be filtered.
    raster_path (str): Path to the raster file.

    Returns:
    gpd.GeoDataFrame: Filtered GeoDataFrame containing only geometries that intersect with the raster.
    """
    assert isinstance(gdf, gpd.GeoDataFrame), "Input gdf must be a GeoDataFrame."

    with rio.open(raster_path) as src:
        if gdf.crs != src.crs:
            gdf = gdf.to_crs(src.crs)

        coordinates = [(geom.x, geom.y) for geom in gdf.geometry]
        gdf['raster_val'] = [val[0] for val in src.sample(coordinates)]

        if src.nodata is not None:
            filtered_gdf = gdf[gdf["raster_val"] != src.nodata]
        else:
            # If no formal nodata is set, drop points that fell completely outside 
            # the raster extent
            filtered_gdf = gdf[gdf["raster_val"].notna()]
    return filtered_gdf


def get_daily_means(df, bands, weights=None, error=None):
    """
    Get daily means for given bands from a DataFrame.

    Parameters:
    df (pd.DataFrame): Input DataFrame with a DateTime index.
    bands (list): The bands for which to calculate daily means.
    weights (str, optional): The column name of the weights to be used in the calculation of the daily means.
    error (str, optional): The column name of the error values to be used in the calculation of the daily means.

    Returns:
    pd.DataFrame: DataFrame with daily means for the specified band.
    """
    assert isinstance(df, pd.DataFrame), "Input df must be a pandas DataFrame."
    assert isinstance(df.index, pd.DatetimeIndex), "Input DataFrame must have a DateTime index."
    assert bands != [], "Bands list cannot be empty."

    if error is not None:
        if weights is None:
            weights = 'error_weight'
            df[weights] = 1 / (df[error] ** 2)
        else:
            df[weights] = df[weights] / (df[error] ** 2)

    # Group by date and calculate the mean for the specified bands
    if weights is not None:
        weighted_bands = df[bands].multiply(df[weights], axis=0)
        weighted_bands[weights] = df[weights]
        daily_means = weighted_bands.groupby(df.index.date).sum()
        daily_means[bands] = daily_means[bands].div(daily_means[weights], axis=0)
    else:
        daily_means = df.groupby(df.index.date)[bands].mean()
    daily_means.index = pd.to_datetime(daily_means.index)
    
    return daily_means[bands]

def running_mean(df, N, bands):
    """Compute the centered running mean of given bands of a DataFrame.

    Parameters:
    df (pd.DataFrame): Input DataFrame with a DateTime index.
    N (int): Window size in number of days (centered around each date index).
    bands (list): The bands for which to compute the running mean.

    Returns:
    array-like: Array of the same length as x containing the centered running mean.
    """
    assert N > 0, "Window size N must be a positive integer."
    assert isinstance(df, pd.DataFrame), "Input x must be a pandas DataFrame."
    assert isinstance(df.index, pd.DatetimeIndex), "Input DataFrame must have a DateTime index."
    assert bands != [], "Bands list cannot be empty."

    means = (df[bands].rolling(window=f'{N}D', center=True, 
                             min_periods=1, closed='both').mean())
    means = means.groupby(means.index.date).mean()
    means.index = pd.to_datetime(means.index)
    return means

def gaussian_running_mean(df, N, bands, std=None):
    """Compute the centered running mean of given bands of a DataFrame.

    Parameters:
    df (pd.DataFrame): Input DataFrame with a DateTime index.
    N (int): Window size in number of days (centered around each date index).
    bands (list): The bands for which to compute the running mean.

    Returns:
    array-like: Array of the same length as x containing the centered running mean.
    """
    assert (N > 0) and (N % 2 == 1), "Window size N must be a odd positive integer."
    assert isinstance(df, pd.DataFrame), "Input x must be a pandas DataFrame."
    assert isinstance(df.index, pd.DatetimeIndex), "Input DataFrame must have a DateTime index."
    assert bands != [], "Bands list cannot be empty."

    date_range = pd.date_range(df.index.min(), df.index.max(), freq='D')
    df = df.reindex(date_range)
    if std is None:
        std = N / 4
    print(N)
    return df[bands].rolling(window=N, center=True, 
                             min_periods=1, closed='both',
                             win_type='gaussian').mean(std=std)

def mean_phenology(df, bands):
    """Compute the mean phenology of given bands of a DataFrame.

    Parameters:
    df (pd.DataFrame): Input DataFrame with a DateTime index.
    bands (list): The bands for which to compute the mean phenology.

    Returns:
    pd.DataFrame: DataFrame containing the mean phenology for the specified bands.
    """
    assert isinstance(df, pd.DataFrame), "Input x must be a pandas DataFrame."
    assert isinstance(df.index, pd.DatetimeIndex), "Input DataFrame must have a DateTime index."
    assert bands != [], "Bands list cannot be empty."

    return df[bands].groupby(df.index.dayofyear).mean()

def add_years_to_doy_index(df, start, end):
    """Creates a DateTime index from a DataFrame with a day-of-year index by 
       repeating values for each year between start and end.

    Parameters:
    df (pd.DataFrame): Input DataFrame with a day-of-year index.
    start (int): The starting year.
    end (int): The ending year.

    Returns:
    pd.DataFrame: DataFrame with DateTime index covering the specified years.
    """
    assert isinstance(df, pd.DataFrame), "Input x must be a pandas DataFrame."
    assert (df.index.min() >= 1) and (df.index.max() <= 366), "Input DataFrame index must be day-of-year (1-366)."
    assert start <= end, "Start year must be less than or equal to end year."

    year_dfs = []
    for year in range(start, end + 1):
        year_df = df.copy()
        dt = pd.to_datetime(f'{year}-01-01')
        year_df['date'] = year_df.index.to_series().apply(lambda x: dt + pd.Timedelta(f'{x}D'))
        year_df = year_df.set_index('date')
        year_dfs.append(year_df)
    all_year_df = pd.concat(year_dfs)
    return all_year_df

def create_sif_analysis_df(df, bands, start_year, end_year, N, weights=None,
                           error=None,gaussian=False, std=None):
    """Create a DataFrame for SIF analysis with daily means, running mean, 
       and mean phenology.

    Parameters:
    df (pd.DataFrame): Input DataFrame with a DateTime index.
    bands (list): The bands to include in the analysis DataFrame.
    start_year (int): The starting year for the analysis.
    end_year (int): The ending year for the analysis.
    N (int): Window size for the running mean in days.
    gaussian (bool): Whether to use a Gaussian window for the running mean.
    std (float, optional): The standard deviation for the Gaussian window.
    mask_path (str): Path to the mask raster file.
    weights (str, optional): The column name of the weights to be used in the calculation of the daily means.

    Returns:
    pd.DataFrame: DataFrame containing daily means, running mean, and mean phenology for the specified bands.
    """
    assert isinstance(df, pd.DataFrame), "Input x must be a pandas DataFrame."
    assert isinstance(df.index, pd.DatetimeIndex), "Input DataFrame must have a DateTime index."
    assert bands != [], "Bands list cannot be empty."
    assert start_year <= end_year, "Start year must be less than or equal to end year."
    assert N > 0, "Window size N must be a positive integer."

    daily_means = get_daily_means(df, bands, weights=weights, error=error)
    if gaussian:
        running_means = gaussian_running_mean(daily_means, N, bands, std=std)
    else:
        running_means = running_mean(daily_means, N, bands)
    mean_pheno = mean_phenology(running_means, bands)
    mean_pheno = add_years_to_doy_index(mean_pheno, 
                                        start_year, 
                                        end_year)
    # Remove doy 366 of non-leap years to avoid duplicates in the index
    mean_pheno = mean_pheno[~mean_pheno.index.duplicated(keep='last')]

    analysis_df = pd.concat([daily_means, running_means, mean_pheno], 
                            axis=1)
    analysis_df.columns = ([f'daily_mean_{band}' for band in bands]
                           + [f'running_mean_{band}' for band in bands]
                           + [f'mean_pheno_{band}' for band in bands])

    return analysis_df


def get_sif_anomalies(df, bands, weights=None, error=None):
    """Compute SIF anomalies for given bands of a DataFrame.

    Parameters:
    df (pd.DataFrame): Input DataFrame with a DateTime index.
    bands (list): The bands for which to compute the anomalies.
    weights (str, optional): The column name of the weights to be used in the calculation of the daily means.
    error (str, optional): The column name of the error values to be used in the calculation of the daily means.

    Returns:
    pd.DataFrame: DataFrame containing the anomalies for the specified bands.
    """
    assert isinstance(df, pd.DataFrame), "Input x must be a pandas DataFrame."
    assert isinstance(df.index, pd.DatetimeIndex), "Input DataFrame must have a DateTime index."
    assert bands != [], "Bands list cannot be empty."

    daily_means = get_daily_means(df, bands, weights=weights, error=error)
    doy = daily_means.index.dayofyear

    # Duplicate daily means to get a circular effect for mean phenology
    prev_daily = (daily_means.copy()
                    .set_index(pd.to_datetime(2019*1000 + doy, format='%Y%j')))
    curr_daily = (daily_means.copy()
                    .set_index(pd.to_datetime(2020*1000 + doy, format='%Y%j')))
    post_daily = (daily_means.copy()
                    .set_index(pd.to_datetime(2021*1000 + doy, format='%Y%j')))
    doy_sif = pd.concat([prev_daily, curr_daily, post_daily])
    doy_sif = doy_sif.sort_index()

    smooth_sif = running_mean(doy_sif, 7, 'SIF_743_TOAnorm')

    prev_30 = (smooth_sif.shift(-30, freq='D').rolling("30D", min_periods=1, closed='left').mean()
                .reset_index()
                .drop_duplicates())
    prev_30 = prev_30.set_index('index')
    next_30 = (smooth_sif.rolling("30D", min_periods=1, closed='right').mean()
                .reset_index()
                .drop_duplicates())
    next_30 = next_30.set_index('index')
    ratios = prev_30.join(next_30, lsuffix='_prev', rsuffix='_next', how='inner')
    ratios['SIF_743_TOAnorm_ratio'] = ratios['SIF_743_TOAnorm_next'] - ratios['SIF_743_TOAnorm_prev']

    focus_ratios = ratios[ratios.index.year == 2020]
    SoS = focus_ratios['SIF_743_TOAnorm_ratio'].idxmin().day_of_year
    EoS = focus_ratios['SIF_743_TOAnorm_ratio'].idxmax().day_of_year

    daily_means['doy'] = doy
    gs_sif = daily_means[daily_means['doy'].between(SoS+15, EoS-15)].copy()

    print(len(gs_sif))
    slope, intercept, _, _ = theilslopes(gs_sif['SIF_743_TOAnorm'].values, gs_sif['doy'].values, method='joint')

    gs_sif['SIF_743_TOAnorm_pred'] = slope * gs_sif['doy'] + intercept
    gs_sif['SIF_743_TOAnorm_residual'] = gs_sif['SIF_743_TOAnorm'] - gs_sif['SIF_743_TOAnorm_pred']

    anoms = {}
    for year in gs_sif.index.year.unique():
        gs_sif_year = gs_sif[gs_sif.index.year == year]
        gs_anom = gs_sif_year[gs_sif_year['doy'].between(160, 210)]
        anoms[year] = gs_anom['SIF_743_TOAnorm_residual'].mean()

    return anoms