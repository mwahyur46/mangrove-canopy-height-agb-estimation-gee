/* ============================================================================
 *  Mangrove Canopy Height + Above-Ground Biomass + Carbon Stock Estimation
 *  West Kalimantan Mangrove Coast
 *
 *  Two-stage RF regression:
 *    Stage 1 — Canopy Height (CH): GEDI L2A RH98 as target
 *              Predictors: S2 + indices (NDVI/NDWI/MNDWI/NDMI/NDRE/SAVI/EVI)
 *              Output: wall-to-wall CH map (masked to GMW mangrove)
 *
 *    Stage 2 — AGB: GEDI L4A AGBD as target
 *              Predictors: S2 + indices (NDVI/NDWI/MNDWI/NDMI/NDRE/SAVI/EVI) + CH map (Stage 1)
 *              Output: wall-to-wall AGB map (masked to GMW mangrove)
 *
 *    Stage 3 — Carbon stock: derived as AGB x 0.451 (IPCC Wetlands Supplement 2013)
 *              Output: wall-to-wall carbon stock map (Mg C/ha)
 *
 *  Author  : Muhammad Wahyu Ramadhan
 *  GitHub  : github.com/mwahyur46
 *  LinkedIn: linkedin.com/in/mwahyur
 * ============================================================================
 *  Required assets (imported via the Code Editor "Imports" panel):
 *    aoi  — ee.Geometry covering West Kalimantan Mangrove Coast
 *
 *  References:
 *    Baloloy, A. B., et al. (2020). Development and application of a new
 *    mangrove vegetation index (MVI) for rapid and accurate mangrove mapping.
 *    ISPRS Journal of Photogrammetry and Remote Sensing, 166, 95-117.
 *    https://doi.org/10.1016/j.isprsjprs.2020.06.001
 *
 *    Duncanson, L., et al. (2022). Aboveground biomass density models for
 *    NASA's Global Ecosystem Dynamics Investigation (GEDI) lidar mission.
 *    Remote Sensing of Environment, 270, 112845.
 *    https://doi.org/10.1016/j.rse.2021.112845
 *
 *    Gupta, K., et al. (2018). An index for discrimination of mangroves from
 *    non-mangroves using LANDSAT 8 OLI imagery. MethodsX, 5, 1129-1139.
 *    https://doi.org/10.1016/j.mex.2018.09.011
 *
 *    Potapov, P., et al. (2020). Mapping global forest canopy height through
 *    integration of GEDI and Landsat data. Remote Sensing of Environment,
 *    253, 112165. https://doi.org/10.1016/j.rse.2020.112165
 *
 *    IPCC (2014). 2013 Supplement to the 2006 IPCC Guidelines for National
 *    Greenhouse Gas Inventories: Wetlands. Chapter 4 (Coastal Wetlands).
 *    Hiraishi, T., Krug, T., Tanabe, K., Srivastava, N., Baasansuren, J.,
 *    Fukuda, M. & Troxler, T.G. (eds). IPCC, Switzerland.
 *    Tier 1 default carbon fraction for mangrove AGB: 0.451 (45.1%).
 *    https://www.ipcc-nggip.iges.or.jp/public/wetlands/
 * ========================================================================== */

// ============================================================================
// 1. CONFIG
// ============================================================================
var S2_START         = '2025-01-01';
var S2_END           = '2025-12-31';
var GEDI_START       = '2019-01-01';
var GEDI_END         = '2025-06-01';
var CLOUD_PCT        = 20;
var SPLIT            = 0.7;
var N_TREES          = 500;
var GMW_ASSET        = 'projects/earthengine-legacy/assets/projects/sat-io/open-datasets/GMW/extent/gmw_v3_2020_vec';
var CARBON_FRACTION  = 0.451;  // IPCC (2014) Wetlands Supplement, Chapter 4 -- mangrove Tier 1 default

Map.centerObject(aoi, 9.5);

// ============================================================================
// 2. SENTINEL-2 — cloud mask + median composite
// ============================================================================
function maskS2clouds(image) {
  var scl  = image.select('SCL');
  var mask = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  return image.updateMask(mask).divide(10000)
              .copyProperties(image, ['system:time_start']);
}

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(S2_START, S2_END)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_PCT))
  .map(maskS2clouds);

print('Sentinel-2 collection (filtered):', s2);

var s2_bands  = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B11', 'B12'];
var s2_median = s2.median().select(s2_bands).clip(aoi);

Map.addLayer(s2_median, {bands: ['B4','B3','B2'], min: 0, max: 0.25, gamma: 0.9},
             'S2 median 2025 (true color)');
Map.addLayer(s2_median, {bands: ['B8','B4','B3'], min: 0, max: 0.4, gamma: 0.9},
             'S2 median 2025 (false color NIR)', false);
Map.addLayer(s2_median, {bands: ['B8A','B11','B4'], min: 0, max: 0.4, gamma: 0.9},
             'S2 median 2025 (false color Red Edge)', false);

// ============================================================================
// 3. SPECTRAL INDICES
// ============================================================================
var ndvi  = s2_median.normalizedDifference(['B8',  'B4' ]).rename('NDVI');
var ndwi  = s2_median.normalizedDifference(['B3',  'B8' ]).rename('NDWI');
var mndwi = s2_median.normalizedDifference(['B3',  'B11']).rename('MNDWI');
var ndmi  = s2_median.normalizedDifference(['B8',  'B11']).rename('NDMI');



// NDRE: Normalized Difference Red Edge
// (B8A - B5) / (B8A + B5) — sensitive to canopy chlorophyll and nitrogen content.
var ndre = s2_median.normalizedDifference(['B8A', 'B5']).rename('NDRE');

// SAVI: Soil Adjusted Vegetation Index
// More robust than NDVI against soil/mudflat background beneath mangrove canopy.
var savi = s2_median.expression(
  '1.5 * (NIR - RED) / (NIR + RED + 0.5)', {
    NIR: s2_median.select('B8'),
    RED: s2_median.select('B4')
  }).rename('SAVI');

// EVI: Enhanced Vegetation Index
// Reduces NDVI saturation in dense mangrove canopy.
var evi = s2_median.expression(
  '2.5 * (NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1)', {
    NIR : s2_median.select('B8'),
    RED : s2_median.select('B4'),
    BLUE: s2_median.select('B2')
  }).rename('EVI');

var indices     = ee.Image.cat([ndvi, ndwi, mndwi, ndmi, ndre, savi, evi]);
var index_bands = ['NDVI', 'NDWI', 'MNDWI', 'NDMI', 'NDRE', 'SAVI', 'EVI'];

Map.addLayer(ndvi, {min: -0.2, max: 0.8, palette: ['white', 'darkgreen']},
             'NDVI', false);

// ============================================================================
// 4. MANGROVE MASK — Global Mangrove Watch v3 (2020)
// ============================================================================
var gmw      = ee.FeatureCollection(GMW_ASSET).filterBounds(aoi);
var gmw_mask = ee.Image(1).clip(gmw).unmask(0).eq(1);

Map.addLayer(gmw_mask.selfMask(), {palette: ['#1a9641']},
             'GMW mangrove mask 2020', false);



// ============================================================================
// 5. BASE FEATURE STACK (shared by both models)
// ============================================================================
// S2 + spectral indices — wall-to-wall, no mask.
// Stack is NOT masked here to prevent null-drop during sampling.
// Slight misalignment between GEDI footprints and the GMW raster mask
// can cause sampleRegions to discard all points if predictors are null.
var base_bands = s2_bands.concat(index_bands);

var base_stack = s2_median
  .addBands(indices)
  .select(base_bands);

print('Base feature stack:', base_stack);

// ============================================================================
// 6. STAGE 1 — CANOPY HEIGHT MODEL
// ============================================================================

// --- 8a. GEDI L2A target ---
// RH98: height at which 98% of returned energy is below.
// quality_flag == 1: good shots only.
function filterL2A(image) {
  return image.updateMask(image.select('quality_flag').eq(1));
}

var gedi_l2a = ee.ImageCollection('LARSE/GEDI/GEDI02_A_002_MONTHLY')
  .filterBounds(aoi)
  .filterDate(GEDI_START, GEDI_END)
  .map(filterL2A);

var rh98_target = gedi_l2a.select('rh98').mean().clip(aoi)
                           .updateMask(gmw_mask);

print('GEDI L2A collection size:', gedi_l2a.size());
print('RH98 valid pixel count:', rh98_target.reduceRegion({
  reducer: ee.Reducer.count(), geometry: aoi, scale: 25, maxPixels: 1e9
}));

Map.addLayer(rh98_target,
             {min: 0, max: 40, palette: ['ffffff','006400']},
             'GEDI RH98 (raw footprints)', false);

// --- 8b. Sampling ---
var ch_training_image = base_stack.addBands(rh98_target.rename('rh98'));

var ch_samples = ch_training_image.sample({
  region    : aoi,
  scale     : 25,
  numPixels : 1e6,
  seed      : 42,
  geometries: true,
  dropNulls : true,
  tileScale : 8
}).randomColumn('rand', 42);

print('CH total samples:', ch_samples.size());

var ch_train = ch_samples.filter(ee.Filter.lt('rand',  SPLIT));
var ch_test  = ch_samples.filter(ee.Filter.gte('rand', SPLIT));

print('CH train samples:', ch_train.size());
print('CH test samples: ', ch_test.size());

// --- 8c. RF regressor ---
var rf_ch = ee.Classifier.smileRandomForest(N_TREES)
              .setOutputMode('REGRESSION')
              .train({
                features       : ch_train,
                classProperty  : 'rh98',
                inputProperties: base_bands
              });

// Wall-to-wall canopy height map
var ch_map = base_stack.classify(rf_ch).rename('CH_m').updateMask(gmw_mask);

var ch_layer = ui.Map.Layer(
  ch_map,
  {min: 10, max: 25, palette: ['ffffb2','fecc5c','fd8d3c','f03b20','bd0026']},
  'Canopy height estimate (m)', false, 1
);
Map.layers().add(ch_layer);

print('Stage 1 (CH) RF trained.');

// --- 8d. CH accuracy assessment ---
var ch_test_pred = ch_test.classify(rf_ch, 'rh98_pred');

var ch_obs      = ch_test_pred.aggregate_array('rh98');
var ch_pred     = ch_test_pred.aggregate_array('rh98_pred');
var ch_n        = ch_test_pred.size();
var ch_obs_mean = ch_obs.reduce(ee.Reducer.mean());

var ch_ss_res = ch_obs.zip(ch_pred).map(function(pair) {
  pair = ee.List(pair);
  var diff = ee.Number(pair.get(0)).subtract(ee.Number(pair.get(1)));
  return diff.pow(2);
}).reduce(ee.Reducer.sum());

var ch_ss_tot = ch_obs.map(function(val) {
  return ee.Number(val).subtract(ch_obs_mean).pow(2);
}).reduce(ee.Reducer.sum());

var ch_r2   = ee.Number(1).subtract(ee.Number(ch_ss_res).divide(ch_ss_tot));
var ch_rmse = ee.Number(ch_ss_res).divide(ch_n).sqrt();
var ch_mae  = ch_obs.zip(ch_pred).map(function(pair) {
  pair = ee.List(pair);
  return ee.Number(pair.get(0)).subtract(ee.Number(pair.get(1))).abs();
}).reduce(ee.Reducer.mean());
var ch_bias = ch_obs.zip(ch_pred).map(function(pair) {
  pair = ee.List(pair);
  return ee.Number(pair.get(1)).subtract(ee.Number(pair.get(0)));
}).reduce(ee.Reducer.mean());

print('--- CH Regression Metrics (test set) ---');
print('R2:  ', ch_r2);
print('RMSE:', ch_rmse, 'm');
print('MAE: ', ch_mae,  'm');
print('Bias:', ch_bias, 'm');

var ch_imp = ee.Dictionary(rf_ch.explain().get('importance'));

// ============================================================================
// 7. STAGE 2 — AGB MODEL
// ============================================================================

// --- 9a. GEDI L4A target ---
var gedi_l4a_col = ee.ImageCollection('LARSE/GEDI/GEDI04_A_002_MONTHLY')
  .filterBounds(aoi)
  .filterDate(GEDI_START, GEDI_END)
  .select(['agbd', 'l4_quality_flag']);

var gedi_l4a_img = gedi_l4a_col.mean().clip(aoi);

var l4a_mask = gedi_l4a_col.select('l4_quality_flag').mean().gt(0.3)
  .and(gedi_l4a_img.select('agbd').gt(0))
  .and(gmw_mask);

var agbd_target = gedi_l4a_img.select('agbd').updateMask(l4a_mask);

print('GEDI L4A collection size:', gedi_l4a_col.size());
print('AGBD valid pixel count:', agbd_target.reduceRegion({
  reducer: ee.Reducer.count(), geometry: aoi, scale: 25, maxPixels: 1e9
}));

Map.addLayer(agbd_target,
             {min: 0, max: 300, palette: ['ffffcc','c2e699','78c679','31a354','006837']},
             'GEDI L4A AGBD (raw footprints)', false);

// --- 9b. AGB feature stack = agb_base + CH map (Stage 1 output) ---
// ch_map is wall-to-wall — safe to use as a predictor band.
var agb_bands = base_bands.concat(['CH_m']);

var agb_stack = base_stack.addBands(ch_map);

print('AGB feature stack (base + CH_m):', agb_stack);

// --- 9c. Sampling ---
var agb_training_image = agb_stack.addBands(agbd_target.rename('agbd'));

var agb_samples = agb_training_image.sample({
  region    : aoi,
  scale     : 25,
  numPixels : 1e6,
  seed      : 42,
  geometries: true,
  dropNulls : true,
  tileScale : 8
}).randomColumn('rand', 42);

print('AGB total samples:', agb_samples.size());

var agb_train = agb_samples.filter(ee.Filter.lt('rand',  SPLIT));
var agb_test  = agb_samples.filter(ee.Filter.gte('rand', SPLIT));

print('AGB train samples:', agb_train.size());
print('AGB test samples: ', agb_test.size());

// --- 9d. RF regressor ---
var rf_agb = ee.Classifier.smileRandomForest(N_TREES)
               .setOutputMode('REGRESSION')
               .train({
                 features       : agb_train,
                 classProperty  : 'agbd',
                 inputProperties: agb_bands
               });

var agb_map = agb_stack.classify(rf_agb).rename('AGB_Mgha').updateMask(gmw_mask);

var agb_layer = ui.Map.Layer(
  agb_map,
  {min: 50, max: 200, palette: ['ffffcc','c2e699','78c679','31a354','006837']},
  'AGB estimate (Mg/ha)', false, 1
);
Map.layers().add(agb_layer);

print('Stage 2 (AGB) RF trained.');

// --- 9e. AGB accuracy assessment ---
var agb_test_pred = agb_test.classify(rf_agb, 'agbd_pred');

var agb_obs      = agb_test_pred.aggregate_array('agbd');
var agb_pred_arr = agb_test_pred.aggregate_array('agbd_pred');
var agb_n        = agb_test_pred.size();
var agb_obs_mean = agb_obs.reduce(ee.Reducer.mean());

var agb_ss_res = agb_obs.zip(agb_pred_arr).map(function(pair) {
  pair = ee.List(pair);
  var diff = ee.Number(pair.get(0)).subtract(ee.Number(pair.get(1)));
  return diff.pow(2);
}).reduce(ee.Reducer.sum());

var agb_ss_tot = agb_obs.map(function(val) {
  return ee.Number(val).subtract(agb_obs_mean).pow(2);
}).reduce(ee.Reducer.sum());

var agb_r2   = ee.Number(1).subtract(ee.Number(agb_ss_res).divide(agb_ss_tot));
var agb_rmse = ee.Number(agb_ss_res).divide(agb_n).sqrt();
var agb_mae  = agb_obs.zip(agb_pred_arr).map(function(pair) {
  pair = ee.List(pair);
  return ee.Number(pair.get(0)).subtract(ee.Number(pair.get(1))).abs();
}).reduce(ee.Reducer.mean());
var agb_bias = agb_obs.zip(agb_pred_arr).map(function(pair) {
  pair = ee.List(pair);
  return ee.Number(pair.get(1)).subtract(ee.Number(pair.get(0)));
}).reduce(ee.Reducer.mean());

print('--- AGB Regression Metrics (test set) ---');
print('R2:  ', agb_r2);
print('RMSE:', agb_rmse, 'Mg/ha');
print('MAE: ', agb_mae,  'Mg/ha');
print('Bias:', agb_bias, 'Mg/ha');

var agb_imp = ee.Dictionary(rf_agb.explain().get('importance'));

// ============================================================================
// STAGE 3 — CARBON STOCK (derived from AGB)
// ============================================================================
// Linear conversion: Carbon (Mg C/ha) = AGB (Mg/ha) x 0.451
//
// Carbon fraction reference:
//   IPCC (2014). 2013 Supplement to the 2006 IPCC Guidelines for National
//   Greenhouse Gas Inventories: Wetlands. Chapter 4 (Coastal Wetlands).
//   Tier 1 default carbon fraction for mangrove AGB = 0.451 (45.1%).
//   This is the gold-standard reference used by Global Mangrove Watch and
//   national greenhouse gas inventories for blue carbon accounting.
var carbon_stock = agb_map.multiply(CARBON_FRACTION)
                           .rename('Carbon_MgC_ha')
                           .updateMask(gmw_mask);

var carbon_layer = ui.Map.Layer(
  carbon_stock,
  {min: 25, max: 90, palette: ['f7f4f9','d4b9da','c994c7','df65b0','980043']},
  'Carbon stock (Mg C/ha)', true, 1
);
Map.layers().add(carbon_layer);

print('Stage 3 (Carbon stock) derived.');

// Carbon stock distribution stats
var carbon_stats = carbon_stock.reduceRegion({
  reducer  : ee.Reducer.mean()
              .combine({reducer2: ee.Reducer.stdDev(),        sharedInputs: true})
              .combine({reducer2: ee.Reducer.minMax(),        sharedInputs: true})
              .combine({reducer2: ee.Reducer.percentile([50]), sharedInputs: true})
              .combine({reducer2: ee.Reducer.sum(),           sharedInputs: true}),
  geometry : aoi,
  scale    : 25,
  maxPixels: 1e10,
  tileScale: 8
});

print('Carbon stock distribution stats:', carbon_stats);

// Total carbon: pixel sum x pixel area (25m x 25m = 625 m2 = 0.0625 ha)
var carbon_total_Mg = ee.Number(carbon_stats.get('Carbon_MgC_ha_sum'))
                       .multiply(0.0625);

print('Total carbon stock (Mg C):', carbon_total_Mg);

// ============================================================================
// 8. CONSOLE CHARTS
// ============================================================================
print(ui.Chart.feature.byFeature({
  features: ch_test_pred, xProperty: 'rh98', yProperties: ['rh98_pred']
}).setChartType('ScatterChart').setOptions({
  title: 'Observed vs Predicted Canopy Height (test set)',
  hAxis: {title: 'Observed RH98 (m)'}, vAxis: {title: 'Predicted RH98 (m)'},
  legend: {position: 'none'}, pointSize: 4, colors: ['#1a73e8'],
  trendlines: {0: {type: 'linear', color: '#e53935', lineWidth: 1.5, opacity: 0.8}}
}));

print(ui.Chart.feature.byFeature({
  features: agb_test_pred, xProperty: 'agbd', yProperties: ['agbd_pred']
}).setChartType('ScatterChart').setOptions({
  title: 'Observed vs Predicted AGBD (test set)',
  hAxis: {title: 'Observed AGBD (Mg/ha)'}, vAxis: {title: 'Predicted AGBD (Mg/ha)'},
  legend: {position: 'none'}, pointSize: 4, colors: ['#1a73e8'],
  trendlines: {0: {type: 'linear', color: '#e53935', lineWidth: 1.5, opacity: 0.8}}
}));

// ============================================================================
// 9. LEFT SIDEBAR
// ============================================================================
var leftPanel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style : {width: '300px', padding: '12px', backgroundColor: 'white'}
});

leftPanel.add(ui.Label('Mangrove CH + AGB + Carbon Stock',
              {fontWeight: 'bold', fontSize: '16px', margin: '0 0 2px 0'}));
leftPanel.add(ui.Label('West Kalimantan Mangrove Coast — 2025',
              {fontSize: '12px', color: '#444', margin: '0 0 2px 0'}));
leftPanel.add(ui.Label('S2 + spectral indices + GEDI L2A/L4A | Two-stage RF Regressor',
              {fontSize: '11px', color: '#666', margin: '0 0 6px 0'}));
leftPanel.add(ui.Label('Muhammad Wahyu Ramadhan',
              {fontSize: '11px', margin: '0'}));
leftPanel.add(ui.Label('github.com/mwahyur46',
              {fontSize: '11px', color: '#1a73e8', margin: '0 0 4px 0'}));
leftPanel.add(ui.Label('', {height: '1px', backgroundColor: '#ccc',
              margin: '8px 0', stretch: 'horizontal'}));

leftPanel.add(ui.Label('Layers',
              {fontWeight: 'bold', fontSize: '13px', margin: '0 0 4px 0'}));

for (var i = 0; i < Map.layers().length(); i++) {
  (function(idx) {
    var layer = Map.layers().get(idx);
    var cb = ui.Checkbox({
      label: layer.getName(), value: layer.getShown(),
      style: {fontSize: '11px', margin: '2px 0'}
    });
    cb.onChange(function(checked) { layer.setShown(checked); });
    leftPanel.add(cb);
  })(i);
}

leftPanel.add(ui.Label('', {height: '1px', backgroundColor: '#ccc',
              margin: '8px 0', stretch: 'horizontal'}));

// Carbon stock legend -- continuous horizontal gradient (purple sequential)
leftPanel.add(ui.Label('Carbon Stock (Mg C/ha)',
              {fontWeight: 'bold', fontSize: '12px', margin: '0 0 4px 0'}));
leftPanel.add(ui.Thumbnail({
  image : ee.Image.pixelLonLat().select('longitude').unitScale(-180, 180)
             .visualize({min: 0, max: 1,
                         palette: ['f7f4f9','d4b9da','c994c7','df65b0','980043']}),
  params: {bbox: [-180, -1, 180, 1], dimensions: '256x16'},
  style : {stretch: 'horizontal', height: '16px', margin: '0', padding: '0'}
}));
leftPanel.add(ui.Panel([
  ui.Label('25',  {fontSize: '10px', margin: '2px 0'}),
  ui.Label('45',  {fontSize: '10px', margin: '2px 0', stretch: 'horizontal', textAlign: 'center'}),
  ui.Label('65',  {fontSize: '10px', margin: '2px 0', stretch: 'horizontal', textAlign: 'center'}),
  ui.Label('90+', {fontSize: '10px', margin: '2px 0'})
], ui.Panel.Layout.flow('horizontal'), {stretch: 'horizontal'}));

leftPanel.add(ui.Label('', {height: '1px', backgroundColor: '#ccc',
              margin: '8px 0', stretch: 'horizontal'}));

// AGB legend
leftPanel.add(ui.Label('AGB (Mg/ha)',
              {fontWeight: 'bold', fontSize: '12px', margin: '0 0 4px 0'}));
leftPanel.add(ui.Thumbnail({
  image : ee.Image.pixelLonLat().select('longitude').unitScale(-180, 180)
             .visualize({min: 0, max: 1,
                         palette: ['ffffcc','c2e699','78c679','31a354','006837']}),
  params: {bbox: [-180, -1, 180, 1], dimensions: '256x16'},
  style : {stretch: 'horizontal', height: '16px', margin: '0', padding: '0'}
}));
leftPanel.add(ui.Panel([
  ui.Label('50',   {fontSize: '10px', margin: '2px 0'}),
  ui.Label('100',  {fontSize: '10px', margin: '2px 0', stretch: 'horizontal', textAlign: 'center'}),
  ui.Label('150',  {fontSize: '10px', margin: '2px 0', stretch: 'horizontal', textAlign: 'center'}),
  ui.Label('200+', {fontSize: '10px', margin: '2px 0'})
], ui.Panel.Layout.flow('horizontal'), {stretch: 'horizontal'}));

leftPanel.add(ui.Label('', {height: '1px', backgroundColor: '#ccc',
              margin: '8px 0', stretch: 'horizontal'}));

// Canopy height legend
leftPanel.add(ui.Label('Canopy Height (m)',
              {fontWeight: 'bold', fontSize: '12px', margin: '0 0 4px 0'}));
leftPanel.add(ui.Thumbnail({
  image : ee.Image.pixelLonLat().select('longitude').unitScale(-180, 180)
             .visualize({min: 0, max: 1,
                         palette: ['ffffb2','fecc5c','fd8d3c','f03b20','bd0026']}),
  params: {bbox: [-180, -1, 180, 1], dimensions: '256x16'},
  style : {stretch: 'horizontal', height: '16px', margin: '0', padding: '0'}
}));
leftPanel.add(ui.Panel([
  ui.Label('10',  {fontSize: '10px', margin: '2px 0'}),
  ui.Label('15',  {fontSize: '10px', margin: '2px 0', stretch: 'horizontal', textAlign: 'center'}),
  ui.Label('20',  {fontSize: '10px', margin: '2px 0', stretch: 'horizontal', textAlign: 'center'}),
  ui.Label('25+', {fontSize: '10px', margin: '2px 0'})
], ui.Panel.Layout.flow('horizontal'), {stretch: 'horizontal'}));

leftPanel.add(ui.Label('', {height: '1px', backgroundColor: '#ccc',
              margin: '8px 0', stretch: 'horizontal'}));

// Opacity sliders
leftPanel.add(ui.Label('Carbon Stock Layer Opacity',
              {fontWeight: 'bold', fontSize: '12px', margin: '0 0 2px 0'}));
leftPanel.add(ui.Slider({
  min: 0, max: 1, value: 1, step: 0.05,
  onChange: function(v) { carbon_layer.setOpacity(v); },
  style: {stretch: 'horizontal'}
}));

leftPanel.add(ui.Label('AGB Layer Opacity',
              {fontWeight: 'bold', fontSize: '12px', margin: '0 0 2px 0'}));
leftPanel.add(ui.Slider({
  min: 0, max: 1, value: 1, step: 0.05,
  onChange: function(v) { agb_layer.setOpacity(v); },
  style: {stretch: 'horizontal'}
}));

leftPanel.add(ui.Label('Canopy Height Layer Opacity',
              {fontWeight: 'bold', fontSize: '12px', margin: '6px 0 2px 0'}));
leftPanel.add(ui.Slider({
  min: 0, max: 1, value: 1, step: 0.05,
  onChange: function(v) { ch_layer.setOpacity(v); },
  style: {stretch: 'horizontal'}
}));

leftPanel.add(ui.Label('', {height: '1px', backgroundColor: '#ccc',
              margin: '10px 0', stretch: 'horizontal'}));

// Pixel inspector
leftPanel.add(ui.Label('Pixel Inspector',
              {fontWeight: 'bold', fontSize: '13px', margin: '0 0 4px 0'}));

var inspectorContent = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'), style: {margin: '0'}
});
inspectorContent.add(ui.Label('Click anywhere on the map to inspect.',
              {fontSize: '11px', color: '#888'}));
leftPanel.add(inspectorContent);

Map.style().set('cursor', 'crosshair');
Map.onClick(function(coords) {
  inspectorContent.clear();
  inspectorContent.add(ui.Label(
    'Lon: ' + coords.lon.toFixed(5) + '  Lat: ' + coords.lat.toFixed(5),
    {fontSize: '10px', color: '#555', margin: '0 0 4px 0'}
  ));
  inspectorContent.add(ui.Label('Sampling...', {fontSize: '11px', color: '#888'}));

  var point     = ee.Geometry.Point([coords.lon, coords.lat]);
  var sampleImg = agb_stack
    .addBands(ch_map)
    .addBands(agb_map)
    .addBands(carbon_stock);

  sampleImg.reduceRegion({
    reducer: ee.Reducer.first(), geometry: point, scale: 25
  }).evaluate(function(vals) {
    inspectorContent.clear();
    inspectorContent.add(ui.Label(
      'Lon: ' + coords.lon.toFixed(5) + '  Lat: ' + coords.lat.toFixed(5),
      {fontSize: '10px', color: '#555', margin: '0 0 4px 0'}
    ));
    if (!vals || vals['AGB_Mgha'] === null || vals['AGB_Mgha'] === undefined) {
      inspectorContent.add(ui.Label('No mangrove data at this location.',
                    {fontSize: '11px', color: '#c00'}));
      return;
    }
    inspectorContent.add(ui.Label(
      'Carbon stock: ' + (vals['Carbon_MgC_ha'] !== null ? vals['Carbon_MgC_ha'].toFixed(2) : 'n/a') + ' Mg C/ha',
      {fontSize: '12px', fontWeight: 'bold', color: '#980043', margin: '0 0 2px 0'}
    ));
    inspectorContent.add(ui.Label(
      'AGB: ' + vals['AGB_Mgha'].toFixed(2) + ' Mg/ha',
      {fontSize: '12px', fontWeight: 'bold', color: '#31a354', margin: '0 0 2px 0'}
    ));
    inspectorContent.add(ui.Label(
      'Canopy height: ' + (vals['CH_m'] !== null ? vals['CH_m'].toFixed(2) : 'n/a') + ' m',
      {fontSize: '12px', fontWeight: 'bold', color: '#1d91c0', margin: '0 0 4px 0'}
    ));
    agb_bands.forEach(function(b) {
      var v = vals[b];
      var display = (v === null || v === undefined) ? 'n/a' :
                    (Math.abs(v) < 0.001 && v !== 0 ? v.toExponential(2) : v.toFixed(4));
      inspectorContent.add(ui.Label(b + ': ' + display,
                    {fontSize: '10px', margin: '1px 0'}));
    });
  });
});

ui.root.insert(0, leftPanel);

// ============================================================================
// 10. RIGHT SIDEBAR — model performance (AGB top, CH bottom)
// ============================================================================
var rightPanel = ui.Panel({
  layout: ui.Panel.Layout.flow('vertical'),
  style : {width: '300px', padding: '12px', backgroundColor: 'white'}
});

rightPanel.add(ui.Label('Model Performance',
              {fontWeight: 'bold', fontSize: '17px', margin: '0 0 2px 0'}));
rightPanel.add(ui.Label('Two-stage RF Regressor | 30% held-out split.',
              {fontSize: '11px', color: '#666', margin: '0 0 4px 0'}));
rightPanel.add(ui.Label('', {height: '1px', backgroundColor: '#ccc',
              margin: '6px 0', stretch: 'horizontal'}));

function r2Color(v) {
  if (v === null || v === undefined || isNaN(v)) return '#888';
  if (v < 0.6) return '#c0392b';
  if (v < 0.8) return '#e67e22';
  return '#27ae60';
}

function addRegressionBlock(parent, title, r2_ee, rmse_ee, mae_ee, bias_ee,
                             unit, test_fc, x_prop, y_prop,
                             x_label, y_label, imp_dict,
                             n_trees, n_train, n_test,
                             gedi_period, s2_period, ref_note) {
  parent.add(ui.Label(title,
              {fontWeight: 'bold', fontSize: '14px', margin: '4px 0 4px 0',
               color: '#222'}));

  // Metrics
  var r2Label   = ui.Label('R2: computing...',   {fontSize: '12px', margin: '2px 0'});
  var rmseLabel = ui.Label('RMSE: computing...', {fontSize: '12px', margin: '2px 0'});
  var maeLabel  = ui.Label('MAE: computing...',  {fontSize: '12px', margin: '2px 0'});
  var biasLabel = ui.Label('Bias: computing...', {fontSize: '12px', margin: '2px 0'});

  parent.add(r2Label);
  parent.add(rmseLabel);
  parent.add(maeLabel);
  parent.add(biasLabel);

  r2_ee.evaluate(function(v) {
    r2Label.setValue('R2: ' + v.toFixed(4));
    r2Label.style().set('color', r2Color(v));
  });
  rmse_ee.evaluate(function(v) { rmseLabel.setValue('RMSE: ' + v.toFixed(3) + ' ' + unit); });
  mae_ee.evaluate(function(v)  { maeLabel.setValue( 'MAE:  ' + v.toFixed(3) + ' ' + unit); });
  bias_ee.evaluate(function(v) { biasLabel.setValue('Bias: ' + v.toFixed(3) + ' ' + unit); });

  // Training info
  parent.add(ui.Label('Training Info',
              {fontWeight: 'bold', fontSize: '12px', margin: '6px 0 2px 0'}));
  parent.add(ui.Label('n Trees    : ' + n_trees,   {fontSize: '11px', margin: '1px 0'}));
  parent.add(ui.Label('Train split: ' + (SPLIT * 100) + '%', {fontSize: '11px', margin: '1px 0'}));
  var trainLabel = ui.Label('Train n    : computing...', {fontSize: '11px', margin: '1px 0'});
  var testLabel  = ui.Label('Test n     : computing...', {fontSize: '11px', margin: '1px 0'});
  parent.add(trainLabel);
  parent.add(testLabel);
  n_train.evaluate(function(v) { trainLabel.setValue('Train n    : ' + v); });
  n_test.evaluate(function(v)  { testLabel.setValue( 'Test n     : ' + v); });
  parent.add(ui.Label('GEDI period: ' + gedi_period, {fontSize: '11px', margin: '1px 0'}));
  parent.add(ui.Label('S2 period  : ' + s2_period,   {fontSize: '11px', margin: '1px 0'}));
  parent.add(ui.Label(ref_note,
              {fontSize: '10px', color: '#666', margin: '1px 0 4px 0'}));

  // Scatter plot
  parent.add(ui.Label('Observed vs Predicted',
              {fontWeight: 'bold', fontSize: '12px', margin: '6px 0 2px 0'}));
  parent.add(ui.Chart.feature.byFeature({
    features: test_fc, xProperty: x_prop, yProperties: [y_prop]
  }).setChartType('ScatterChart').setOptions({
    title    : '',
    hAxis    : {title: x_label, textStyle: {fontSize: 9}},
    vAxis    : {title: y_label, textStyle: {fontSize: 9}},
    legend   : {position: 'none'},
    pointSize: 3,
    colors   : ['#1a73e8'],
    trendlines: {0: {type: 'linear', color: '#e53935', lineWidth: 1.5,
                     opacity: 0.8, visibleInLegend: false}},
    width    : 270, height: 180,
    chartArea: {left: 45, top: 10, right: 10, bottom: 45}
  }));

  // Feature importance
  parent.add(ui.Label('Feature Importance',
              {fontWeight: 'bold', fontSize: '12px', margin: '6px 0 2px 0'}));
  parent.add(ui.Chart.array.values({
    array: imp_dict.values(), axis: 0, xLabels: imp_dict.keys()
  }).setChartType('ColumnChart').setOptions({
    title    : '',
    legend   : {position: 'none'},
    width    : 270, height: 160,
    hAxis    : {textStyle: {fontSize: 8}, slantedText: true, slantedTextAngle: 45},
    vAxis    : {textStyle: {fontSize: 9}},
    chartArea: {left: 40, top: 10, right: 10, bottom: 55}
  }));
}

// ============================================================
// Carbon Stock block (top) -- derived, not modeled
// ============================================================
rightPanel.add(ui.Label('Carbon Stock Estimation (Stage 3)',
              {fontWeight: 'bold', fontSize: '14px', margin: '4px 0 4px 0', color: '#222'}));
rightPanel.add(ui.Label('Derived: AGB x 0.451 (IPCC carbon fraction)',
              {fontSize: '11px', color: '#666', margin: '0 0 2px 0'}));

var cMean   = ui.Label('Mean      : computing...', {fontSize: '12px', margin: '2px 0'});
var cMedian = ui.Label('Median    : computing...', {fontSize: '12px', margin: '2px 0'});
var cStd    = ui.Label('Std Dev   : computing...', {fontSize: '12px', margin: '2px 0'});
var cMin    = ui.Label('Min       : computing...', {fontSize: '12px', margin: '2px 0'});
var cMax    = ui.Label('Max       : computing...', {fontSize: '12px', margin: '2px 0'});
var cTotal  = ui.Label('Total C   : computing...', {fontSize: '12px', margin: '2px 0', fontWeight: 'bold'});

rightPanel.add(cMean);
rightPanel.add(cMedian);
rightPanel.add(cStd);
rightPanel.add(cMin);
rightPanel.add(cMax);
rightPanel.add(cTotal);

carbon_stats.evaluate(function(s) {
  if (!s) return;
  cMean.setValue(  'Mean      : ' + s['Carbon_MgC_ha_mean'].toFixed(2)    + ' Mg C/ha');
  cMedian.setValue('Median    : ' + s['Carbon_MgC_ha_p50'].toFixed(2)     + ' Mg C/ha');
  cStd.setValue(   'Std Dev   : ' + s['Carbon_MgC_ha_stdDev'].toFixed(2)  + ' Mg C/ha');
  cMin.setValue(   'Min       : ' + s['Carbon_MgC_ha_min'].toFixed(2)     + ' Mg C/ha');
  cMax.setValue(   'Max       : ' + s['Carbon_MgC_ha_max'].toFixed(2)     + ' Mg C/ha');
});
carbon_total_Mg.evaluate(function(v) {
  cTotal.setValue('Total C   : ' + (v / 1000).toFixed(2) + ' kt C');
});

rightPanel.add(ui.Label('Ref: IPCC (2014) Wetlands Supplement, Chapter 4; carbon fraction = 0.451',
              {fontSize: '10px', color: '#666', margin: '4px 0'}));

rightPanel.add(ui.Label('', {height: '1px', backgroundColor: '#ccc',
              margin: '8px 0', stretch: 'horizontal'}));

// AGB block (atas)
addRegressionBlock(
  rightPanel,
  'AGB Estimation (Stage 2)',
  agb_r2, agb_rmse, agb_mae, agb_bias,
  'Mg/ha',
  agb_test_pred, 'agbd', 'agbd_pred',
  'Observed AGBD (Mg/ha)', 'Predicted AGBD (Mg/ha)',
  agb_imp,
  N_TREES, agb_train.size(), agb_test.size(),
  '2019-2025', '2025',
  'Ref: GEDI L4A (Duncanson et al. 2022)'
);

rightPanel.add(ui.Label('', {height: '1px', backgroundColor: '#ccc',
              margin: '8px 0', stretch: 'horizontal'}));

// CH block (bottom)
addRegressionBlock(
  rightPanel,
  'Canopy Height (Stage 1)',
  ch_r2, ch_rmse, ch_mae, ch_bias,
  'm',
  ch_test_pred, 'rh98', 'rh98_pred',
  'Observed RH98 (m)', 'Predicted RH98 (m)',
  ch_imp,
  N_TREES, ch_train.size(), ch_test.size(),
  '2019-2025', '2025',
  'Ref: GEDI L2A RH98 (quality_flag == 1); CH mapping adapted from Potapov et al. (2020)'
);

ui.root.add(rightPanel);

// ============================================================================
// 11. EXPORTS
// ============================================================================
var EXPORT_FOLDER = 'GEE_agb_ch_west_kalimantan';

Export.image.toDrive({
  image         : agb_map.toFloat(),
  description   : 'agb_map_west_kalimantan_2025',
  fileNamePrefix: 'agb_map_west_kalimantan_2025',
  folder        : EXPORT_FOLDER,
  region        : aoi, scale: 25, maxPixels: 1e10, fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image         : ch_map.toFloat(),
  description   : 'canopy_height_map_west_kalimantan_2025',
  fileNamePrefix: 'canopy_height_map_west_kalimantan_2025',
  folder        : EXPORT_FOLDER,
  region        : aoi, scale: 25, maxPixels: 1e10, fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image         : agb_stack.toFloat(),
  description   : 'feature_stack_agb_west_kalimantan_2025',
  fileNamePrefix: 'feature_stack_agb_west_kalimantan_2025',
  folder        : EXPORT_FOLDER,
  region        : aoi, scale: 25, maxPixels: 1e10, fileFormat: 'GeoTIFF'
});

Export.table.toDrive({
  collection    : agb_test_pred,
  description   : 'agb_test_predictions_west_kalimantan',
  fileNamePrefix: 'agb_test_predictions_west_kalimantan',
  folder        : EXPORT_FOLDER, fileFormat: 'CSV'
});

Export.table.toDrive({
  collection    : ch_test_pred,
  description   : 'ch_test_predictions_west_kalimantan',
  fileNamePrefix: 'ch_test_predictions_west_kalimantan',
  folder        : EXPORT_FOLDER, fileFormat: 'CSV'
});

Export.image.toDrive({
  image         : carbon_stock.toFloat(),
  description   : 'carbon_stock_map_west_kalimantan_2025',
  fileNamePrefix: 'carbon_stock_map_west_kalimantan_2025',
  folder        : EXPORT_FOLDER,
  region        : aoi, scale: 25, maxPixels: 1e10, fileFormat: 'GeoTIFF'
});