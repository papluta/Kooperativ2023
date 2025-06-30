var roi = ee.Geometry.Rectangle([9.52282041, 51.5340748, 10.21825130, 51.89208507]);

//Terra and Aqua LST Daily Global 1km
var modis = ee.ImageCollection("MODIS/006/MOD11A1")
                  .filterBounds(roi)
                  .filter(ee.Filter.date('2000-01-01', '2023-06-30'))
  .select(['LST_Day_1km', 'LST_Night_1km', 'QC_Day', 'QC_Night']);


// Function to filter and scale both Day and Night LST
var processLST = function(image) {
  var day_raw = image.select('LST_Day_1km');
  var night_raw = image.select('LST_Night_1km');
  var qc_day = image.select('QC_Day');
  var qc_night = image.select('QC_Night');

  // QA mask: keep only pixels where bits 0–1 == 00 (best quality)
  var qa_mask_day = qc_day.bitwiseAnd(3).eq(0);
  var qa_mask_night = qc_night.bitwiseAnd(3).eq(0);


  // Combined masks
  var mask_day = qa_mask_day;
  var mask_night = qa_mask_night;

  // Scale and convert to Celsius
  var lst_day = day_raw.multiply(0.02).subtract(273.15).updateMask(mask_day);
  var lst_night = night_raw.multiply(0.02).subtract(273.15).updateMask(mask_night);

  // Combine: average day and night
  //   var mean_lst = lst_day.add(lst_night).divide(2)
  //       .copyProperties(image, ['system:time_start']);

  // return only day data
  return lst_day;
};

// Apply processing function
var lst_combined = modis.map(processLST);

// Compute mean over time
var lstMean = lst_combined.mean().clip(roi);

// Display
Map.centerObject(roi, 8);
Map.addLayer(lstMean, 
  {min: 0, max: 35, palette: ['blue', 'green', 'yellow', 'orange', 'red']}, 
  'Mean LST 2000–2023');
Map.addLayer(roi, {color: 'black'}, 'ROI');

// Export to Google Drive (optional)
Export.image.toDrive({
  image: lstMean,
  description: 'Kooperativ_MODIS_LST_Mean_2000_2023',
  folder: '',
  region: roi,
  scale: 1000,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});


print(ui.Chart.image.histogram(lstMean, roi, 1000));
