# Mangrove Canopy Height, Above-Ground Biomass, and Carbon Stock Estimation (GEE)

![Platform](https://img.shields.io/badge/platform-Google%20Earth%20Engine-4285F4.svg)
![Language](https://img.shields.io/badge/language-JavaScript-yellow.svg)
![Sensors](https://img.shields.io/badge/sensors-Sentinel--2%20%7C%20GEDI%20L2A%20%7C%20GEDI%20L4A-informational.svg)
![Method](https://img.shields.io/badge/method-Two--stage%20Random%20Forest-green.svg)
![Status](https://img.shields.io/badge/status-completed-green.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

A reproducible **two-stage Random Forest regression** workflow for wall-to-wall mapping of mangrove **canopy height (CH)**, **above-ground biomass density (AGBD)**, and **carbon stock** along the **West Kalimantan Mangrove Coast**, implemented entirely in **Google Earth Engine**. Stage 1 estimates canopy height from GEDI L2A RH98 footprints as training labels; Stage 2 uses the resulting wall-to-wall CH map as an additional predictor for AGBD estimation against GEDI L4A ground truth. Carbon stock is derived from AGB via the IPCC (2014) Wetlands Supplement default carbon fraction for mangrove (0.451). Mangrove extent is constrained using the **Global Mangrove Watch v3 (2020)** mask. The script ships with an interactive dual-sidebar UI (layer toggles, gradient legends, opacity sliders, pixel inspector, and per-model performance panels).

<p align="center">
  <img src="images/map_overview.png" alt="Full map view with carbon stock estimate and dual sidebars" width="95%">
</p>

---

## Background

Accurate estimation of mangrove biomass, canopy structure, and carbon stock is critical for blue carbon accounting, ecosystem monitoring, and conservation planning. Conventional field-based approaches are limited in spatial coverage, while single-stage spectral regression from optical imagery often underperforms due to the spectral saturation of dense canopy. The integration of GEDI lidar with multispectral imagery for canopy height mapping has been demonstrated at global scale by Potapov et al. (2020), who combined GEDI and Landsat to produce a wall-to-wall forest height product. This repository adapts that concept to the mangrove domain using Sentinel-2, implementing a two-stage approach in which GEDI-derived canopy height serves as a structural intermediate predictor for AGBD estimation, substantially improving accuracy over a direct spectral-to-biomass model. Carbon stock is subsequently derived from the AGB map using the Tier 1 default carbon fraction from the IPCC (2014) Wetlands Supplement, Chapter 4 (Coastal Wetlands).

GEDI L4A AGBD values are derived from the model of Duncanson et al. (2022), which converts GEDI RH metrics to above-ground biomass density in Mg/ha. The Global Mangrove Watch v3 extent layer (Bunting et al. 2022) constrains all analysis to confirmed mangrove pixels.

---

## Repository Structure

```
mangrove-canopy-height-agb-estimation-gee/
|-- agb_canopy_height_west_kalimantan.js   # two-stage CH + AGB + carbon stock estimation
|-- images/
|   |-- map_overview.png
|   |-- carbon_stock_map.png
|   |-- agb_map.png
|   |-- canopy_height_map.png
|   |-- false_color_mangrove.png
|   |-- gedi_footprints.png
|   |-- right_panel_carbon.png
|   |-- right_panel_agb.png
|   |-- right_panel_ch.png
|   |-- pixel_inspector_carbon.png
|   |-- pixel_inspector.png
|   |-- scatterplot_agb.png
|   |-- scatterplot_ch.png
|   |-- fimportance_agb.png
|   |-- fimportance_ch.png
|   |-- full_map.png
|   |-- full_right_panel_agb.png
|   |-- full_right_panel_ch.png
|-- README.md
|-- LICENSE
```

---

## Method

The workflow runs entirely server-side in Google Earth Engine and produces three wall-to-wall raster outputs masked to GMW mangrove extent.

| Step | Description |
|------|-------------|
| 1 | Build 2025 annual **median composite** from Sentinel-2 SR Harmonized (SCL-based cloud mask). |
| 2 | Derive **7 spectral indices** (NDVI, NDWI, MNDWI, NDMI, NDRE, SAVI, EVI). |
| 3 | Apply **Global Mangrove Watch v3 (2020)** as binary spatial mask. |
| 4 | Filter **GEDI L2A** (2019-2025, quality_flag == 1); compute mean RH98 per pixel as CH target. |
| 5 | Train **Stage 1 RF regressor** (500 trees, 70/30 split) to predict RH98 from S2 feature stack. |
| 6 | Apply Stage 1 model to produce **wall-to-wall canopy height map**. |
| 7 | Filter **GEDI L4A** (2019-2025, l4_quality_flag > 0.3); compute mean AGBD per pixel. |
| 8 | Train **Stage 2 RF regressor** on S2 feature stack + CH map as additional predictor. |
| 9 | Apply Stage 2 model to produce **wall-to-wall AGB map**. |
| 10 | Derive **carbon stock map**: AGB x 0.451 (IPCC 2014 Wetlands Supplement, Chapter 4). |
| 11 | Evaluate both RF models with R2, RMSE, MAE, and Bias on held-out test set. |
| 12 | Export rasters and test prediction tables to Google Drive. |

### Two-Stage Architecture

```
S2 bands + indices
        |
   [Stage 1 RF]  <-- GEDI L2A RH98 (training labels)
        |
  CH map (wall-to-wall)
        |
S2 bands + indices + CH_m
        |
   [Stage 2 RF]  <-- GEDI L4A AGBD (training labels)
        |
  AGB map (wall-to-wall)
        |
  AGB x 0.451
        |
  Carbon stock map (wall-to-wall)
```

### Spectral Indices

| Index | Formula | Purpose |
|-------|---------|---------|
| NDVI  | (B8 - B4) / (B8 + B4) | Vegetation vigor |
| NDWI  | (B3 - B8) / (B3 + B8) | Open water |
| MNDWI | (B3 - B11) / (B3 + B11) | Water vs built-up separation |
| NDMI  | (B8 - B11) / (B8 + B11) | Canopy moisture |
| NDRE  | (B8A - B5) / (B8A + B5) | Canopy chlorophyll and nitrogen |
| SAVI  | 1.5 * (B8 - B4) / (B8 + B4 + 0.5) | Vegetation with soil adjustment |
| EVI   | 2.5 * (B8 - B4) / (B8 + 6*B4 - 7.5*B2 + 1) | Dense canopy (reduces NDVI saturation) |

### GEDI Products

| Product | Band | Use |
|---------|------|-----|
| LARSE/GEDI/GEDI02_A_002_MONTHLY | rh98 | Canopy height training labels (Stage 1) |
| LARSE/GEDI/GEDI04_A_002_MONTHLY | agbd | AGBD training labels (Stage 2) |

---

## How to Run

1. Open the script in the [Google Earth Engine Code Editor](https://code.earthengine.google.com/1492430de8952238d71746f32119d9e9).
2. Copy `agb_canopy_height_west_kalimantan.js` into a new GEE script.
3. Import the following asset through the **Imports** panel:
   - `aoi` -- `ee.Geometry` covering the target mangrove area
4. Click **Run**. Stage 1 and Stage 2 models train sequentially. The interactive sidebars populate after server-side computation completes (typically 30 to 60 seconds depending on AOI size).
5. Submit the **Export** tasks from the Tasks panel to push outputs to Google Drive.

> **Note on reproducibility.** Google Earth Engine's Random Forest implementation (`smileRandomForest`) does not expose a fixed random seed for tree building. Model performance metrics may vary slightly across runs due to stochastic variation in the RF ensemble and server-side tile partitioning during sampling. Results reported here represent a representative run. Larger deviations may occur during periods of high GEE server load. For publication-grade reproducibility, exporting the sample CSV and retraining in Python (scikit-learn with `random_state=42`) is recommended.

---

## Interactive UI

The script renders a docked dual-sidebar layout eliminating the need to navigate the GEE Console for metrics.

### Left Panel: Layer Controls, Legends, and Pixel Inspector

Layer checkboxes toggle each map layer. Continuous gradient legends are shown for Carbon Stock (Mg C/ha, purple), AGB (Mg/ha, green), and Canopy Height (m, yellow-red). Three independent opacity sliders allow transparency control of each estimate layer. The pixel inspector populates on map click with carbon stock, AGB, canopy height, and all predictor band values.

<p align="center">
  <img src="images/pixel_inspector_carbon.png" alt="Pixel inspector showing carbon stock, AGB, canopy height, and band values" width="40%">
</p>

### Right Panel: Model Performance

Three blocks are shown in the right panel -- Carbon Stock (Stage 3) at the top, AGB (Stage 2) in the middle, and CH (Stage 1) at the bottom. The carbon stock block reports distributional statistics (mean, median, std dev, min, max, total carbon in kt C). Each RF model block shows R2, RMSE, MAE, Bias, training info, scatter plot, and feature importance chart.

<p align="center">
  <img src="images/right_panel_carbon.png" alt="Right panel: Carbon Stock Estimation block" width="32%">
  <img src="images/right_panel_agb.png" alt="Right panel: AGB model performance" width="32%">
  <img src="images/right_panel_ch.png" alt="Right panel: CH model performance" width="32%">
</p>

---

## Data Requirements

| Source | Product | Period | Resolution | Role |
|--------|---------|--------|------------|------|
| Copernicus | Sentinel-2 SR Harmonized | 2025-01-01 to 2025-12-31 | 10 m | Spectral predictors |
| NASA / LARSE | GEDI L2A Monthly | 2019-01-01 to 2025-06-01 | ~25 m footprint | CH training labels |
| NASA / LARSE | GEDI L4A Monthly | 2019-01-01 to 2025-06-01 | ~25 m footprint | AGBD training labels |
| Bunting et al. | Global Mangrove Watch v3 2020 | 2020 | 25 m | Spatial mask |

All datasets are accessed directly via the GEE data catalog. No local downloads are required.

AOI: West Kalimantan Mangrove Coast (south of Pontianak, Kalimantan Barat, Indonesia).

---

## Reproducibility

| Order | Script | Description |
|-------|--------|-------------|
| 1 | `agb_canopy_height_west_kalimantan.js` | Two-stage CH + AGB estimation + carbon stock derivation |

Run the script in GEE Code Editor with the `aoi` asset imported.

---

## Results

### Output Maps

<p align="center">
  <img src="images/carbon_stock_map.png" alt="Wall-to-wall carbon stock estimate (Mg C/ha)" width="95%">
</p>

<p align="center">
  <img src="images/agb_map.png" alt="Wall-to-wall AGB estimate (Mg/ha)" width="48%">
  <img src="images/canopy_height_map.png" alt="Wall-to-wall canopy height estimate (m)" width="48%">
</p>

False color composite (B8A/B11/B4) highlighting mangrove extent:

<p align="center">
  <img src="images/false_color_mangrove.png" alt="False color composite B8A/B11/B4 highlighting mangrove" width="95%">
</p>

GEDI footprint distribution within the AOI (diagonal orbit pattern):

<p align="center">
  <img src="images/gedi_footprints.png" alt="GEDI L2A and L4A footprint distribution in AOI" width="95%">
</p>

### Validation Metrics

Results from a representative run. Due to stochastic variation in GEE's RF implementation, R2 values may vary across runs -- see the note under How to Run.

| Model | R2 | RMSE | MAE | Bias |
|-------|----|------|-----|------|
| Stage 1: Canopy Height (RH98) | 0.3273 | 6.300 m | 4.651 m | -0.129 m |
| Stage 2: AGB (AGBD) | 0.7279 | 38.790 Mg/ha | 23.956 Mg/ha | 0.960 Mg/ha |

### Carbon Stock Statistics (AOI-wide)

Derived from AGB map using IPCC (2014) carbon fraction = 0.451. No separate model validation applies -- uncertainty propagates directly from the AGB model.

| Statistic | Value |
|-----------|-------|
| Mean | 47.33 Mg C/ha |
| Median | 41.00 Mg C/ha |
| Std Dev | 20.03 Mg C/ha |
| Min | 3.90 Mg C/ha |
| Max | 286.20 Mg C/ha |
| **Total Carbon** | **5,588.79 kt C** |

### Scatter Plots: Observed vs Predicted

<p align="center">
  <img src="images/scatterplot_ch.png" alt="Observed vs predicted canopy height (RH98)" width="48%">
  <img src="images/scatterplot_agb.png" alt="Observed vs predicted AGBD" width="48%">
</p>

### Feature Importance

<p align="center">
  <img src="images/fimportance_ch.png" alt="Feature importance: Stage 1 canopy height model" width="48%">
  <img src="images/fimportance_agb.png" alt="Feature importance: Stage 2 AGB model" width="48%">
</p>

For the CH model, B6 (red-edge, 740 nm) is the dominant predictor, consistent with its sensitivity to canopy chlorophyll content and vertical structure. For the AGB model, CH_m is the single most important predictor by a large margin, validating the two-stage design.

### Caveats

- **Model variability.** GEE's `smileRandomForest` does not support a fixed random seed for tree construction. Results may vary across runs. For publication-grade reproducibility, exporting the sample CSV and retraining in Python (scikit-learn with `random_state=42`) is recommended.
- **GEDI footprint sparsity.** GEDI coverage follows a diagonal orbital pattern. Areas between orbit tracks are estimated by the RF model, not directly observed by GEDI. Prediction uncertainty is higher in under-sampled zones.
- **Temporal mismatch.** Sentinel-2 composite uses 2025 imagery; GEDI labels span 2019-2025. Mangrove structural changes over this period are not accounted for.
- **Carbon stock uncertainty.** The 0.451 carbon fraction is a Tier 1 global default. Species-specific or region-specific values may differ. Uncertainty in the AGB estimate propagates linearly to the carbon stock estimate.

### Enhanced Feature Set

Adding Sentinel-1 SAR (VV/VH), SRTM slope, CMRI, and MVI to the feature stack improves model performance. The left panel subtitle reflects the expanded sensor stack.

| Model | R2 | RMSE | MAE | Bias |
|-------|----|------|-----|------|
| Stage 1: Canopy Height (RH98) | 0.6119 | 5.024 m | 3.267 m | -0.059 m |
| Stage 2: AGB (AGBD) | 0.7113 | 41.415 Mg/ha | 25.041 Mg/ha | 1.490 Mg/ha |

<p align="center">
  <img src="images/full_map.png" alt="Full map view: enhanced version with AGB estimate" width="95%">
</p>

<p align="center">
  <img src="images/full_right_panel_agb.png" alt="Enhanced version: AGB model performance panel" width="45%">
  <img src="images/full_right_panel_ch.png" alt="Enhanced version: CH model performance panel" width="45%">
</p>

In the enhanced version feature importance, VV ranks among the top predictors for CH alongside B6, confirming that SAR backscatter captures canopy structural information beyond what optical bands provide. For AGB, CH_m remains the dominant predictor, with VH and slope contributing meaningful secondary signal.

---

## Citation

If this work supports your research or project:

> Ramadhan, M. W. (2026). *Mangrove Canopy Height, Above-Ground Biomass, and Carbon Stock Estimation Using Two-Stage Random Forest and GEDI in Google Earth Engine*. GitHub repository. https://github.com/mwahyur46/mangrove-canopy-height-agb-estimation-gee

---

## References

- Baloloy, A. B., Blanco, A. C., Ana, R. R. C. S., & Nadaoka, K. (2020). Development and application of a new mangrove vegetation index (MVI) for rapid and accurate mangrove mapping. ISPRS Journal of Photogrammetry and Remote Sensing, 166, 95-117. https://doi.org/10.1016/j.isprsjprs.2020.06.001
- Bunting, P., Rosenqvist, A., Hilarides, L., Lucas, R. M., Thomas, N., Tadono, T., ... & Rebelo, L. M. (2022). Global mangrove extent change 1996-2020: Global Mangrove Watch Version 3.0. Remote Sensing, 14(15), 3657. https://doi.org/10.3390/rs14153657
- Duncanson, L., Kellner, J. R., Armston, J., Dubayah, R., Disney, M., Healey, S. P., ... & Hancock, S. (2022). Aboveground biomass density models for NASA's Global Ecosystem Dynamics Investigation (GEDI) lidar mission. Remote Sensing of Environment, 270, 112845. https://doi.org/10.1016/j.rse.2021.112845
- Gupta, K., Mukhopadhyay, A., Giri, S., Chanda, A., Majumdar, S. D., Samanta, S., ... & Hazra, S. (2018). An index for discrimination of mangroves from non-mangroves using LANDSAT 8 OLI imagery. MethodsX, 5, 1129-1139. https://doi.org/10.1016/j.mex.2018.09.011
- IPCC (2014). 2013 Supplement to the 2006 IPCC Guidelines for National Greenhouse Gas Inventories: Wetlands. Chapter 4 (Coastal Wetlands). Hiraishi, T., Krug, T., Tanabe, K., Srivastava, N., Baasansuren, J., Fukuda, M. & Troxler, T.G. (eds). IPCC, Switzerland. https://www.ipcc-nggip.iges.or.jp/public/wetlands/
- Main-Knorn, M., Pflug, B., Louis, J., Debaecker, V., Muller-Wilm, U., & Gascon, F. (2017). Sen2Cor for Sentinel-2. Image and Signal Processing for Remote Sensing XXIII, 10427, 37-48. https://doi.org/10.1117/12.2278218
- Potapov, P., Li, X., Hernandez-Serna, A., Tyukavina, A., Hansen, M. C., Kommareddy, A., Pickens, A., Turubanova, S., Tang, H., Silva, C. E., Armston, J., Dubayah, R., Blair, J. B., & Hofton, M. (2020). Mapping global forest canopy height through integration of GEDI and Landsat data. Remote Sensing of Environment, 253, 112165. https://doi.org/10.1016/j.rse.2020.112165

---

## Acknowledgements

This repository is a personal portfolio project developed to demonstrate applied geospatial data science methodology. Academic background: Master of Remote Sensing, Faculty of Geography, Universitas Gadjah Mada.

Datasets are provided by NASA (GEDI), the European Space Agency (Sentinel-2), and the Global Mangrove Watch consortium. Processing was conducted on the Google Earth Engine platform.

---

## Contact

- Muhammad Wahyu Ramadhan
- GitHub: github.com/mwahyur46
- LinkedIn: linkedin.com/in/mwahyur
- Email: mwahyur46@gmail.com
