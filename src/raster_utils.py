"""
Raster I/O and wall-to-wall inference utilities.

Handles loading the multi-band predictor stack exported from GEE, applying
a trained scikit-learn model pixel-by-pixel to produce a prediction raster,
masking to the mangrove extent, and writing results to GeoTIFF.
"""

import numpy as np
import rasterio
from rasterio.plot import show


def load_stack(raster_path):
    """
    Load a multi-band GeoTIFF predictor stack into a (rows, cols, bands)
    array, along with its rasterio profile and band descriptions.

    Args:
        raster_path (str): Path to the GeoTIFF stack.

    Returns:
        tuple: (stack_array, profile, band_names)
            stack_array (np.ndarray): Shape (height, width, n_bands).
            profile (dict): Rasterio profile for writing outputs.
            band_names (list[str]): Band descriptions, if set in the file.
    """
    with rasterio.open(raster_path) as src:
        data = src.read()  # (bands, height, width)
        profile = src.profile.copy()
        band_names = list(src.descriptions)

    stack_array = np.moveaxis(data, 0, -1)  # (height, width, bands)

    print(f'Loaded stack: {raster_path}')
    print(f'  Shape (H, W, bands): {stack_array.shape}')
    print(f'  Bands: {band_names}')

    return stack_array, profile, band_names


def stack_to_dataframe(stack_array, band_names):
    """
    Flatten a (H, W, bands) raster stack into a 2D feature table for
    prediction, preserving the valid-pixel mask for reconstruction.

    Args:
        stack_array (np.ndarray): Shape (height, width, n_bands).
        band_names (list[str]): Column names matching band order.

    Returns:
        tuple: (X_flat, valid_mask, original_shape)
            X_flat (np.ndarray): Shape (n_valid_pixels, n_bands).
            valid_mask (np.ndarray): Boolean, shape (height, width).
            original_shape (tuple): (height, width).
    """
    height, width, n_bands = stack_array.shape
    flat = stack_array.reshape(-1, n_bands)

    # Valid pixel: no NaN across any band.
    valid_mask_flat = ~np.isnan(flat).any(axis=1)
    valid_mask = valid_mask_flat.reshape(height, width)

    X_flat = flat[valid_mask_flat]

    print(f'Total pixels: {height * width}')
    print(f'Valid pixels: {X_flat.shape[0]}')

    return X_flat, valid_mask, (height, width)


def predict_to_raster(model, X_flat, valid_mask, original_shape, nodata=np.nan):
    """
    Apply a trained model to flattened valid pixels and reconstruct a
    2D prediction raster.

    Args:
        model: Fitted estimator with `.predict()`.
        X_flat (np.ndarray): Shape (n_valid_pixels, n_features).
        valid_mask (np.ndarray): Boolean, shape (height, width).
        original_shape (tuple): (height, width).
        nodata (float): Fill value for invalid pixels.

    Returns:
        np.ndarray: Prediction raster, shape (height, width).
    """
    predictions = model.predict(X_flat)

    output = np.full(original_shape, nodata, dtype='float32')
    output[valid_mask] = predictions

    return output


def apply_mangrove_mask(raster_array, mask_array, nodata=np.nan):
    """
    Mask a prediction raster to the mangrove extent.

    Args:
        raster_array (np.ndarray): Shape (height, width).
        mask_array (np.ndarray): Binary mask, same shape, 1 = keep.
        nodata (float): Fill value outside the mask.

    Returns:
        np.ndarray: Masked raster, shape (height, width).
    """
    masked = raster_array.copy()
    masked[mask_array != 1] = nodata
    return masked


def write_geotiff(array, profile, output_path, nodata=np.nan, dtype='float32'):
    """
    Write a single-band 2D array to GeoTIFF using a reference profile.

    Args:
        array (np.ndarray): Shape (height, width).
        profile (dict): Rasterio profile (from load_stack), used as the
            geospatial reference. Band count and dtype are overridden.
        output_path (str): Output file path.
        nodata (float): NoData value to set in metadata.
        dtype (str): Output data type.
    """
    out_profile = profile.copy()
    out_profile.update(count=1, dtype=dtype, nodata=nodata)

    with rasterio.open(output_path, 'w', **out_profile) as dst:
        dst.write(array.astype(dtype), 1)

    print(f'Written: {output_path}')


def compute_carbon_stock(agb_array, carbon_fraction=0.451):
    """
    Derive carbon stock from an AGB raster via deterministic raster math.

    Carbon (Mg C/ha) = AGB (Mg/ha) x carbon_fraction

    Args:
        agb_array (np.ndarray): AGB raster, shape (height, width).
        carbon_fraction (float): IPCC (2014) Tier 1 mangrove carbon fraction.

    Returns:
        np.ndarray: Carbon stock raster, same shape as input.
    """
    return agb_array * carbon_fraction


def raster_summary_stats(array, pixel_area_ha=None, label='Raster'):
    """
    Print distributional summary statistics for a raster, ignoring NaN.

    Args:
        array (np.ndarray): Raster array.
        pixel_area_ha (float, optional): Pixel area in hectares. If given,
            also prints the total sum scaled to that area.
        label (str): Label for the printed header.

    Returns:
        dict: {'mean', 'median', 'std', 'min', 'max', 'total'(optional)}
    """
    valid = array[~np.isnan(array)]

    stats = {
        'mean': float(np.mean(valid)),
        'median': float(np.median(valid)),
        'std': float(np.std(valid)),
        'min': float(np.min(valid)),
        'max': float(np.max(valid)),
    }

    print(f'--- {label} Distribution Stats ---')
    print(f"Mean  : {stats['mean']:.2f}")
    print(f"Median: {stats['median']:.2f}")
    print(f"StdDev: {stats['std']:.2f}")
    print(f"Min   : {stats['min']:.2f}")
    print(f"Max   : {stats['max']:.2f}")

    if pixel_area_ha is not None:
        total = float(np.sum(valid)) * pixel_area_ha
        stats['total'] = total
        print(f'Total : {total:.2f} (sum x {pixel_area_ha} ha/pixel)')

    return stats
