# The code in this module was generated with the assistance of 
# Anthropic's Claude Sonnet 5 model. Code was then reviewed by 
# Cameron Scholl to ensure it met the requirements of the project.

import numpy as np

def compute_local_semivariogram(array, mask, pixel_width, focal_coords,
                                 lag_width, max_lag=None):
    assert array.shape == mask.shape, "array and mask must have the same shape"
    assert pixel_width > 0, "pixel_width must be positive"
    assert lag_width > 0, "lag_width must be positive"

    n_rows, n_cols = array.shape
    focal_row, focal_col = focal_coords
    assert 0 <= focal_row < n_rows and 0 <= focal_col < n_cols, \
        "focal_coords must fall within the array"
    assert mask[focal_row, focal_col], "the focal pixel must be valid"

    if max_lag is None:
        max_lag = min(focal_row, n_rows - 1 - focal_row,
                       focal_col, n_cols - 1 - focal_col) * pixel_width

    rows, cols = np.indices(array.shape)
    distances = np.hypot(rows - focal_row, cols - focal_col) * pixel_width

    # Valid pixels are not-masked, within maximum lag, and not the focal pixel.
    valid = mask & (distances > 0) & (distances <= max_lag)
    distances = distances[valid]
    squared_diffs = (array[valid] - array[focal_row, focal_col]) ** 2

    if distances.size == 0:
        print("No valid pixels found within the specified lag distance.")
        empty = np.array([])
        return empty, empty, empty

    bin_indices = np.floor_divide(distances, lag_width, dtype=np.int32)

    lag_distances = []
    semivariances = []
    n_pairs = []
    for bin_index in range(max_lag//lag_width + 1):
        in_bin = bin_indices == bin_index
        lag_distances.append(distances[in_bin].mean())
        semivariances.append(0.5 * squared_diffs[in_bin].mean())
        n_pairs.append(in_bin.sum())

    return np.array(lag_distances), np.array(semivariances), np.array(n_pairs)
