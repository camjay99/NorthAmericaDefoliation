import pandas as pd

def get_daily_means(df, bands):
    """
    Get daily means for given bands from a DataFrame.

    Parameters:
    df (pd.DataFrame): Input DataFrame with a DateTime index.
    bands (list): The bands for which to calculate daily means.

    Returns:
    pd.DataFrame: DataFrame with daily means for the specified band.
    """
    assert isinstance(df, pd.DataFrame), "Input df must be a pandas DataFrame."
    assert isinstance(df.index, pd.DatetimeIndex), "Input DataFrame must have a DateTime index."
    assert bands != [], "Bands list cannot be empty."

    # Group by date and calculate the mean for the specified bands
    daily_means = df.groupby(df.index.date)[bands].mean()
    daily_means.index = pd.to_datetime(daily_means.index)
    
    return daily_means

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

    return df[bands].rolling(window=f'{N}D', center=True, 
                             min_periods=1, closed='both').mean()

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

def create_sif_analysis_df(df, bands, start_year, end_year, N):
    """Create a DataFrame for SIF analysis with daily means, running mean, 
       and mean phenology.

    Parameters:
    df (pd.DataFrame): Input DataFrame with a DateTime index.
    bands (list): The bands to include in the analysis DataFrame.
    start_year (int): The starting year for the analysis.
    end_year (int): The ending year for the analysis.
    N (int): Window size for the running mean in days.

    Returns:
    pd.DataFrame: DataFrame containing daily means, running mean, and mean phenology for the specified bands.
    """
    assert isinstance(df, pd.DataFrame), "Input x must be a pandas DataFrame."
    assert isinstance(df.index, pd.DatetimeIndex), "Input DataFrame must have a DateTime index."
    assert bands != [], "Bands list cannot be empty."
    assert start_year <= end_year, "Start year must be less than or equal to end year."
    assert N > 0, "Window size N must be a positive integer."

    daily_means = get_daily_means(df, bands)
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