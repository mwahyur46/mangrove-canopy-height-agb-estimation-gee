"""
Google Earth Engine utilities for mangrove canopy height, AGB, and carbon
stock estimation.

This module is a Python (earthengine-api) port of the GEE JavaScript
"full" workflow, used to prepare the predictor stack, sample GEDI training
data, and export rasters/tables for local modeling in scikit-learn.

Feature set (full version):
    S2 bands         : B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
    Spectral indices : NDVI, NDWI, MNDWI, NDMI, CMRI, MVI, NDRE, SAVI, EVI
    SAR               : VV, VH (Sentinel-1, speckle-filtered)
    Topography        : slope (SRTM)

References:
    Baloloy, A. B., et al. (2020). Development and application of a new
    mangrove vegetation index (MVI) for rapid and accurate mangrove mapping.
    ISPRS Journal of Photogrammetry and Remote Sensing, 166, 95-117.
    https://doi.org/10.1016/j.isprsjprs.2020.06.001

    Duncanson, L., et al. (2022). Aboveground biomass density models for
    NASA's Global Ecosystem Dynamics Investigation (GEDI) lidar mission.
    Remote Sensing of Environment, 270, 112845.
    https://doi.org/10.1016/j.rse.2021.112845

    Gupta, K., et al. (2018). An index for discrimination of mangroves from
    non-mangroves using LANDSAT 8 OLI imagery. MethodsX, 5, 1129-1139.
    https://doi.org/10.1016/j.mex.2018.09.011

    Potapov, P., et al. (2020). Mapping global forest canopy height through
    integration of GEDI and Landsat data. Remote Sensing of Environment,
    253, 112165. https://doi.org/10.1016/j.rse.2020.112165

    IPCC (2014). 2013 Supplement to the 2006 IPCC Guidelines for National
    Greenhouse Gas Inventories: Wetlands. Chapter 4 (Coastal Wetlands).
    Tier 1 default carbon fraction for mangrove AGB: 0.451 (45.1%).
    https://www.ipcc-nggip.iges.or.jp/public/wetlands/
"""

import ee


# ============================================================================
# CONFIG (defaults, override as needed when calling functions)
# ============================================================================
S2_BANDS = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12']
INDEX_BANDS = ['NDVI', 'NDWI', 'MNDWI', 'NDMI', 'CMRI', 'MVI', 'NDRE', 'SAVI', 'EVI']
SAR_BANDS = ['VV', 'VH']
TOPO_BANDS = ['slope']

GMW_ASSET = (
    'projects/earthengine-legacy/assets/projects/sat-io/'
    'open-datasets/GMW/extent/gmw_v3_2020_vec'
)
CARBON_FRACTION = 0.451  # IPCC (2014) Wetlands Supplement, mangrove Tier 1 default


# ============================================================================
# 1. SENTINEL-2 COMPOSITE
# ============================================================================
def mask_s2_clouds(image):
    """
    Mask clouds and cloud shadows in a Sentinel-2 SR Harmonized image using
    the Scene Classification Layer (SCL).

    SCL values masked out: 3 (cloud shadow), 8 (cloud medium probability),
    9 (cloud high probability), 10 (thin cirrus).

    Args:
        image (ee.Image): Sentinel-2 SR Harmonized image with SCL band.

    Returns:
        ee.Image: Cloud-masked, reflectance-scaled (divided by 10000) image.
    """
    scl = image.select('SCL')
    mask = scl.neq(3).And(scl.neq(8)).And(scl.neq(9)).And(scl.neq(10))
    return image.updateMask(mask).divide(10000).copyProperties(
        image, ['system:time_start']
    )


def get_s2_median(aoi, start_date, end_date, cloud_pct=20, bands=None):
    """
    Build a cloud-masked Sentinel-2 median composite over an AOI.

    Args:
        aoi (ee.Geometry): Area of interest.
        start_date (str): Start date, 'YYYY-MM-DD'.
        end_date (str): End date, 'YYYY-MM-DD'.
        cloud_pct (float): Max CLOUDY_PIXEL_PERCENTAGE filter threshold.
        bands (list[str], optional): Bands to select. Defaults to S2_BANDS.

    Returns:
        ee.Image: Median composite, clipped to aoi, selected bands only.
    """
    bands = bands or S2_BANDS

    s2 = (
        ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
        .filterBounds(aoi)
        .filterDate(start_date, end_date)
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloud_pct))
        .map(mask_s2_clouds)
    )

    return s2.median().select(bands).clip(aoi)


# ============================================================================
# 2. SPECTRAL INDICES
# ============================================================================
def compute_indices(s2_median):
    """
    Compute the full spectral index set used as predictors.

    Indices: NDVI, NDWI, MNDWI, NDMI, CMRI (Gupta et al. 2018),
    MVI (Baloloy et al. 2020), NDRE, SAVI, EVI.

    Args:
        s2_median (ee.Image): Sentinel-2 median composite with standard bands.

    Returns:
        ee.Image: Multi-band image, one band per index, named per INDEX_BANDS.
    """
    ndvi = s2_median.normalizedDifference(['B8', 'B4']).rename('NDVI')
    ndwi = s2_median.normalizedDifference(['B3', 'B8']).rename('NDWI')
    mndwi = s2_median.normalizedDifference(['B3', 'B11']).rename('MNDWI')
    ndmi = s2_median.normalizedDifference(['B8', 'B11']).rename('NDMI')

    # CMRI: Combined Mangrove Recognition Index (Gupta et al. 2018)
    cmri = ndvi.subtract(ndwi).rename('CMRI')

    # MVI: Mangrove Vegetation Index (Baloloy et al. 2020)
    mvi = s2_median.expression(
        '(NIR - GREEN) / (SWIR1 - GREEN)',
        {
            'NIR': s2_median.select('B8'),
            'GREEN': s2_median.select('B3'),
            'SWIR1': s2_median.select('B11'),
        },
    ).rename('MVI')

    # NDRE: sensitive to canopy chlorophyll and nitrogen content.
    ndre = s2_median.normalizedDifference(['B8A', 'B5']).rename('NDRE')

    # SAVI: robust against soil/mudflat background beneath mangrove canopy.
    savi = s2_median.expression(
        '1.5 * (NIR - RED) / (NIR + RED + 0.5)',
        {'NIR': s2_median.select('B8'), 'RED': s2_median.select('B4')},
    ).rename('SAVI')

    # EVI: reduces NDVI saturation in dense mangrove canopy.
    evi = s2_median.expression(
        '2.5 * (NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1)',
        {
            'NIR': s2_median.select('B8'),
            'RED': s2_median.select('B4'),
            'BLUE': s2_median.select('B2'),
        },
    ).rename('EVI')

    return ee.Image.cat([ndvi, ndwi, mndwi, ndmi, cmri, mvi, ndre, savi, evi])


# ============================================================================
# 3. SENTINEL-1 SAR
# ============================================================================
def get_s1_filtered(aoi, start_date, end_date, smoothing_radius=20):
    """
    Build a speckle-filtered Sentinel-1 VV/VH composite.

    Filters: IW mode, VV+VH polarization, descending orbit, 10 m resolution.
    Speckle reduction via focal mean (circular kernel).

    Args:
        aoi (ee.Geometry): Area of interest.
        start_date (str): Start date, 'YYYY-MM-DD'.
        end_date (str): End date, 'YYYY-MM-DD'.
        smoothing_radius (float): Focal mean radius in meters.

    Returns:
        ee.Image: Two-band image (VV, VH), clipped to aoi.
    """
    s1 = (
        ee.ImageCollection('COPERNICUS/S1_GRD')
        .filterBounds(aoi)
        .filterDate(start_date, end_date)
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        .filter(ee.Filter.eq('resolution_meters', 10))
    )

    vv = (
        s1.select('VV').median().clip(aoi)
        .focal_mean(smoothing_radius, 'circle', 'meters').rename('VV')
    )
    vh = (
        s1.select('VH').median().clip(aoi)
        .focal_mean(smoothing_radius, 'circle', 'meters').rename('VH')
    )

    return ee.Image.cat([vv, vh])


# ============================================================================
# 4. TOPOGRAPHY
# ============================================================================
def get_slope(aoi):
    """
    Compute slope from SRTM 30 m elevation.

    Args:
        aoi (ee.Geometry): Area of interest.

    Returns:
        ee.Image: Single-band 'slope' image (degrees), clipped to aoi.
    """
    return (
        ee.Terrain.slope(ee.Image('USGS/SRTMGL1_003').select('elevation'))
        .clip(aoi)
        .rename('slope')
    )


# ============================================================================
# 5. MANGROVE MASK
# ============================================================================
def get_gmw_mask(aoi, gmw_asset=GMW_ASSET):
    """
    Build a binary mangrove mask from Global Mangrove Watch v3 (2020).

    Args:
        aoi (ee.Geometry): Area of interest.
        gmw_asset (str): GMW FeatureCollection asset path.

    Returns:
        ee.Image: Binary mask (1 = mangrove, 0 = non-mangrove).
    """
    gmw = ee.FeatureCollection(gmw_asset).filterBounds(aoi)
    return ee.Image(1).clip(gmw).unmask(0).eq(1)


# ============================================================================
# 6. BASE FEATURE STACK
# ============================================================================
def build_base_stack(aoi, s2_start, s2_end, cloud_pct=20, smoothing_radius=20):
    """
    Build the full predictor stack: S2 + indices + S1 + slope.

    The stack is NOT masked to the mangrove extent here, to avoid dropping
    GEDI footprints during sampling due to slight raster/footprint
    misalignment at the GMW mask boundary.

    Args:
        aoi (ee.Geometry): Area of interest.
        s2_start (str): Sentinel-2 composite start date.
        s2_end (str): Sentinel-2 composite end date.
        cloud_pct (float): Max cloud percentage filter.
        smoothing_radius (float): S1 speckle filter radius (meters).

    Returns:
        tuple[ee.Image, list[str]]: (base_stack, base_bands)
    """
    s2_median = get_s2_median(aoi, s2_start, s2_end, cloud_pct)
    indices = compute_indices(s2_median)
    s1 = get_s1_filtered(aoi, s2_start, s2_end, smoothing_radius)
    slope = get_slope(aoi)

    base_bands = S2_BANDS + INDEX_BANDS + SAR_BANDS + TOPO_BANDS

    base_stack = (
        s2_median.addBands(indices).addBands(s1).addBands(slope).select(base_bands)
    )

    return base_stack, base_bands


# ============================================================================
# 7. GEDI TARGETS
# ============================================================================
def filter_l2a_quality(image):
    """Keep only good-quality GEDI L2A shots (quality_flag == 1)."""
    return image.updateMask(image.select('quality_flag').eq(1))


def get_ch_target(aoi, gedi_start, gedi_end, gmw_mask):
    """
    Build the canopy height target (GEDI L2A RH98), masked to mangrove extent.

    Args:
        aoi (ee.Geometry): Area of interest.
        gedi_start (str): GEDI collection start date.
        gedi_end (str): GEDI collection end date.
        gmw_mask (ee.Image): Binary mangrove mask.

    Returns:
        tuple[ee.Image, ee.ImageCollection]: (rh98_target, gedi_l2a_collection)
    """
    gedi_l2a = (
        ee.ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY')
        .filterBounds(aoi)
        .filterDate(gedi_start, gedi_end)
        .map(filter_l2a_quality)
    )

    rh98_target = gedi_l2a.select('rh98').mean().clip(aoi).updateMask(gmw_mask)

    return rh98_target, gedi_l2a


def get_agb_target(aoi, gedi_start, gedi_end, gmw_mask, quality_threshold=0.3):
    """
    Build the AGB target (GEDI L4A AGBD), masked to mangrove extent.

    Quality filter: mean l4_quality_flag > threshold AND agbd > 0.

    Args:
        aoi (ee.Geometry): Area of interest.
        gedi_start (str): GEDI collection start date.
        gedi_end (str): GEDI collection end date.
        gmw_mask (ee.Image): Binary mangrove mask.
        quality_threshold (float): Minimum mean l4_quality_flag to keep a pixel.

    Returns:
        tuple[ee.Image, ee.ImageCollection]: (agbd_target, gedi_l4a_collection)
    """
    gedi_l4a_col = (
        ee.ImageCollection('LARSE/GEDI/GEDI04_A_002_MONTHLY')
        .filterBounds(aoi)
        .filterDate(gedi_start, gedi_end)
        .select(['agbd', 'l4_quality_flag'])
    )

    gedi_l4a_img = gedi_l4a_col.mean().clip(aoi)

    l4a_mask = (
        gedi_l4a_col.select('l4_quality_flag').mean().gt(quality_threshold)
        .And(gedi_l4a_img.select('agbd').gt(0))
        .And(gmw_mask)
    )

    agbd_target = gedi_l4a_img.select('agbd').updateMask(l4a_mask)

    return agbd_target, gedi_l4a_col


# ============================================================================
# 8. SAMPLING
# ============================================================================
def sample_training_data(
    image_with_target,
    aoi,
    scale=25,
    num_pixels=int(1e6),
    seed=42,
    tile_scale=8,
):
    """
    Sample a predictor+target image at pixel locations within the AOI.

    Mirrors the GEE JS `ee.Image.sample()` call used for both CH and AGB
    training data extraction. Null predictor/target pixels are dropped.

    Args:
        image_with_target (ee.Image): Predictor stack with target band added.
        aoi (ee.Geometry): Area of interest.
        scale (float): Sampling scale in meters.
        num_pixels (int): Approximate number of pixels to sample.
        seed (int): Random seed for sampling.
        tile_scale (int): Tile scale factor to avoid memory/timeout errors.

    Returns:
        ee.FeatureCollection: Sampled points with predictor + target properties.
    """
    return image_with_target.sample(
        region=aoi,
        scale=scale,
        numPixels=num_pixels,
        seed=seed,
        geometries=True,
        dropNulls=True,
        tileScale=tile_scale,
    )


def export_samples_to_drive(
    feature_collection,
    description,
    folder,
    file_format='CSV',
):
    """
    Start a Drive export task for a sampled FeatureCollection.

    Args:
        feature_collection (ee.FeatureCollection): Samples to export.
        description (str): Task description / base filename.
        folder (str): Google Drive folder name.
        file_format (str): Export file format (default 'CSV').

    Returns:
        ee.batch.Task: Started export task. Check status via task.status().
    """
    task = ee.batch.Export.table.toDrive(
        collection=feature_collection,
        description=description,
        fileNamePrefix=description,
        folder=folder,
        fileFormat=file_format,
    )
    task.start()
    return task


def export_image_to_drive(
    image,
    description,
    aoi,
    folder,
    scale=25,
    max_pixels=int(1e10),
):
    """
    Start a Drive export task for an image (e.g. the base feature stack).

    Args:
        image (ee.Image): Image to export.
        description (str): Task description / base filename.
        aoi (ee.Geometry): Export region.
        folder (str): Google Drive folder name.
        scale (float): Export resolution in meters.
        max_pixels (int): Max pixel count allowed for export.

    Returns:
        ee.batch.Task: Started export task. Check status via task.status().
    """
    task = ee.batch.Export.image.toDrive(
        image=image.toFloat(),
        description=description,
        fileNamePrefix=description,
        folder=folder,
        region=aoi,
        scale=scale,
        maxPixels=max_pixels,
        fileFormat='GeoTIFF',
    )
    task.start()
    return task
