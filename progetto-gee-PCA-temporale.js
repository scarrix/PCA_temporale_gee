//------------------------------------------------------------------------------
// 0. DEFINIZIONE ROI E PARAMETRI DI ANALISI
//------------------------------------------------------------------------------
var roi = ee.FeatureCollection(vigneto);  // Sostituisci con la tua ROI vigneto
Map.centerObject(roi, 16);
Map.addLayer(roi, {color: 'green'}, 'Ulivo');

var startYear = 2017;
var endYear   = 2024;

var calculateAndPrintROIStats = function(roi) {
  // Calcola l'area in metri quadrati. Il parametro '1' è l'errore massimo consentito.
  var area_sq_m = roi.geometry().area(1);
  // Converte l'area in ettari (1 ettaro = 10,000 m^2).
  var area_ha = area_sq_m.divide(10000);
  
  // Per contare i pixel, creiamo un'immagine costante e usiamo reduceRegion.
  var pixel_count = ee.Image.constant(1).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: roi.geometry(),
    scale: 10, // Scala coerente con i dati Sentinel (10m).
    maxPixels: 1e13
  }).get('constant'); // Il risultato è in un dizionario con la chiave del nome della banda.

  // Usa .evaluate() per recuperare i risultati dal server in modo asincrono
  // e stamparli nella console del browser.
  ee.Dictionary({area: area_ha, pixels: pixel_count}).evaluate(function(results, error) {
    if (error) {
      print('Errore nel calcolo delle statistiche:', error);
    } else {
      print('Statistiche ROI:',
            'Superficie: ' + results.area.toFixed(2) + ' ettari',
            'Numero di pixel (a 10m): ' + results.pixels);
    }
  });
};

// Esegui la funzione per ciascuna ROI.
calculateAndPrintROIStats(roi);

//------------------------------------------------------------------------------
// 1. PREPARAZIONE DATI SENTINEL-1 CON INDICI
//------------------------------------------------------------------------------

var bands     = ['VV', 'VH', 'RVI', 'RFDI'];

//LEE ANTI-SPECKLE FILTER
function applyLeeFilter(image) {
  var vv = ee.Image(image).select('VV'),
      vh = ee.Image(image).select('VH');

  //Convert dB to linear
  var vvLin = ee.Image(10).pow(vv.divide(10)),
      vhLin = ee.Image(10).pow(vh.divide(10));

  var kernel = ee.Kernel.square({radius: 1});

  function leeOneBand(bandLin, name) {
    var mean = bandLin.reduceNeighborhood(ee.Reducer.mean(), kernel);
    var variance = bandLin.reduceNeighborhood(ee.Reducer.variance(), kernel);
    var cv = variance.sqrt().divide(mean);
    var k = cv.multiply(cv).divide(cv.multiply(cv).add(0.25));
    var filtered = mean.add(k.multiply(bandLin.subtract(mean)));
    // back to dB
    return ee.Image(10).multiply(filtered.log10()).rename(name);
  }

  var vvF = leeOneBand(vvLin, 'VV'),
      vhF = leeOneBand(vhLin, 'VH');

  return image.addBands([vvF, vhF], null, true);
}


function calculateRVI(image) {
  var vh = image.select('VH');
  var vv = image.select('VV');
  
  var rvi = vv.divide(vv.add(vh)).sqrt()
               .multiply(vv.divide(vh)).rename('RVI');
  return image.addBands(rvi);
}

function calculateRFDI(image) {
  var vvLin = ee.Image(10).pow(image.select('VV').divide(10));
  var vhLin = ee.Image(10).pow(image.select('VH').divide(10));
  var rfdi = vvLin.subtract(vhLin)
                  .divide(vvLin.add(vhLin))
                  .rename('RFDI');
  return image.addBands(rfdi);
}

function addAllIndices(image) {
  // Apply speckle filter
  image = applyLeeFilter(image);
  image = calculateRVI(image);
  image = calculateRFDI(image);
  var date = ee.Date(image.get('system:time_start'))
                 .format('YYYY-MM-dd HH:mm:ss');
  return image.set('datetime', date);
}

//Carica Sentinel-1 e aggiungi indici
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(roi)
  .filterDate(startYear + '-01-01', endYear + '-12-31')
  .filterMetadata('transmitterReceiverPolarisation','equals',['VV','VH'])
  .filterMetadata('instrumentMode','equals','IW')
  .filter(ee.Filter.eq('orbitProperties_pass','ASCENDING'))
  .map(addAllIndices);

//------------------------------------------------------------------------------
// 2. SERIE TEMPORALI DEGLI INDICI
//------------------------------------------------------------------------------
print('--- Grafici Serie Storiche ---');
print(
  ui.Chart.image.series({
    imageCollection: s1.select(['VV','VH']),
    region: roi, reducer: ee.Reducer.mean(), scale: 30, xProperty: 'datetime'
  })
  .setChartType('ScatterChart')
  .setOptions({
    title: 'Backscatter (dB) VV & VH',
    series: {0:{label:'VV'},1:{label:'VH'}},
    lineWidth:1, pointSize:3,
    hAxis:{title:'Data'}, vAxis:{title:'Backscatter (dB)'}
  })
);
print(
  ui.Chart.image.series({
    imageCollection: s1.select(['RVI', 'RFDI']),
    region: roi, reducer: ee.Reducer.mean(), scale: 30, xProperty: 'datetime'
  })
  .setChartType('ScatterChart')
  .setOptions({
    title: 'Indici Radar: RVI & RFDI',
    series: {0:{label:'RVI'},1:{label:'RFDI'}},
    lineWidth:1, pointSize:3,
    hAxis:{title:'Data'}
  })
);

//--------------------------------------------------------------------------------
// 3. Calcola medie annuali, per anno, delle varie bande(VV,VH,RVI,RFDI) della ROI
//--------------------------------------------------------------------------------

print('--- Normalizzazione Z-score e PCA ---');

function annualMean(year) {
  var img = s1.filter(ee.Filter.calendarRange(year, year, 'year'))
             .select(bands)
             .mean();
  var stats = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: roi,
    scale: 10,
    maxPixels: 1e9
  });
  return ee.Feature(null, stats).set('year', year);
}
var years = ee.List.sequence(startYear, endYear);
var annualStats = ee.FeatureCollection(years.map(annualMean));

//Estrai valori mediani per ogni banda
var vvValues   = annualStats.aggregate_array('VV');
var vhValues   = annualStats.aggregate_array('VH');
var rviValues  = annualStats.aggregate_array('RVI');
var rfdiValues = annualStats.aggregate_array('RFDI');

print('VV Values:', vvValues);
print('VH Values:', vhValues);
print('RVI Values:', rviValues);
print('RFDI Values:', rfdiValues);

//------------------------------------------------------------------------------
// 3. ANALISI DELLE COMPONENTI PRINCIPALI (PCA)
//------------------------------------------------------------------------------

// 3.1 Normalizzazione (Z-Score)
var vvMean = vvValues.reduce(ee.Reducer.mean());
var vhMean = vhValues.reduce(ee.Reducer.mean());
var rviMean = rviValues.reduce(ee.Reducer.mean());
var rfdiMean = rfdiValues.reduce(ee.Reducer.mean());

var vvStd = vvValues.reduce(ee.Reducer.stdDev());
var vhStd = vhValues.reduce(ee.Reducer.stdDev());
var rviStd = rviValues.reduce(ee.Reducer.stdDev());
var rfdiStd = rfdiValues.reduce(ee.Reducer.stdDev());

print('Means - VV:', vvMean, 'VH:', vhMean, 'RVI:', rviMean, 'RFDI:', rfdiMean);
print('StdDev - VV:', vvStd, 'VH:', vhStd, 'RVI:', rviStd, 'RFDI:', rfdiStd);

var vvNormalized = vvValues.map(function(val) { 
  return ee.Number(val).subtract(vvMean).divide(vvStd); 
});
var vhNormalized = vhValues.map(function(val) { 
  return ee.Number(val).subtract(vhMean).divide(vhStd); 
});
var rviNormalized = rviValues.map(function(val) { 
  return ee.Number(val).subtract(rviMean).divide(rviStd); 
});
var rfdiNormalized = rfdiValues.map(function(val) { 
  return ee.Number(val).subtract(rfdiMean).divide(rfdiStd); 
});

// 3.2 Costruisci matrice normalizzata (già centrata e scalata)
var dataMatrix = ee.Array([vvNormalized, vhNormalized, rviNormalized, rfdiNormalized]).transpose();
print('Normalized data matrix:', dataMatrix);

// Calcola matrice di correlazione (matrice di covarianza)
// Poiché i dati sono già normalizzati, la matrice di covarianza è equivalente alla matrice di correlazione
var n = ee.Number(years.length());
var corrMatrix = dataMatrix.transpose().matrixMultiply(dataMatrix).divide(n);
print('Correlation matrix:', corrMatrix);

// 3.3 Decomposizione autovettori/autovalori
var eigens = corrMatrix.eigen();
var eigenValues = eigens.slice(1, 0, 1).project([0]);
var eigenVectors = eigens.slice(1, 1);
print('Eigenvalues:', eigenValues);
print('Eigenvectors:', eigenVectors);

//Gli autovettori definiscono le “direzioni” nello spazio delle variabili (VV, VH, RVI, RFDI) lungo le quali la varianza è massima.
//Gli autovalori corrispondenti ti dicono quanta varianza (percentuale) è catturata da ciascuna di queste direzioni.

//------------------------------------------------------------------------------
// 4. VISUALIZZAZIONE RISULTATI PCA
//------------------------------------------------------------------------------

print('--- Visualizzazione dei risultati! ---');

// Proiezione sui PC (score) (punteggi/score ottenuti proiettando i tuoi dati standardizzati (gli Z‑score annuali) sugli autovettori)
var pcPerYears = dataMatrix.matrixMultiply(eigenVectors);
print('PC per anno:', pcPerYears); //mostra la matrice completa years x Pc (years = righe, PC = colonne)

// Costruisci FeatureCollection per visualizzazione
var pcFeatures = ee.FeatureCollection(
  ee.List.sequence(0, years.length().subtract(1)).map(function(idx) {
    var year = years.get(idx);
    var pc1 = pcPerYears.slice(0, idx, ee.Number(idx).add(1))
                   .slice(1, 0, 1)
                   .project([0])
                   .get([0]);
    var pc2 = pcPerYears.slice(0, idx, ee.Number(idx).add(1))
                   .slice(1, 1, 2)
                   .project([0])
                   .get([0]);
    return ee.Feature(null, {
      'year': year,
      'PC1': pc1,
      'PC2': pc2
    });
  })
);
// print('PC Features:', pcFeatures);// Mostra una FeatureCollection per ogni anno. Ogni anno ha come proprietà numeroAnno, PC1, PC2.

// Grafico PC1 e PC2 nel tempo (score in funzione del tempo)
var chart = ui.Chart.feature.byFeature(pcFeatures, 'year', ['PC1', 'PC2'])
  .setChartType('LineChart')
  .setOptions({
    title: 'Evoluzione Temporale PC1 e PC2',
    hAxis: {title: 'Anno'},
    vAxis: {title: 'Score PC'},
    lineWidth: 2,
    pointSize: 4,
    series: {
      0: {label: 'PC1', color: 'blue'},
      1: {label: 'PC2', color: 'red'}
    },
    legend: {position: 'right'}
  });
print(chart);

// Calcola varianza spiegata
var totalVariance = eigenValues.reduce(ee.Reducer.sum(), [0]).get([0]);
var varianceExplained = eigenValues.divide(ee.Number(totalVariance));
print('Total Variance:', totalVariance);
print('Variance explained:', varianceExplained);


// Scree Plot - Varianza spiegata da ciascuna componente principale
var eigenList = eigenValues.toList();
var varianceList = varianceExplained.toList();


// Crea FeatureCollection per lo scree plot
var screeFeatures = ee.FeatureCollection([
  ee.Feature(null, {
    'Component': 'PC1',
    'ComponentNum': 1,
    'Eigenvalue': eigenList.get(0),
    'VarianceExplained': varianceList.get(0)
  }),
  ee.Feature(null, {
    'Component': 'PC2',
    'ComponentNum': 2,
    'Eigenvalue': eigenList.get(1),
    'VarianceExplained': varianceList.get(1)
  }),
  ee.Feature(null, {
    'Component': 'PC3',
    'ComponentNum': 3,
    'Eigenvalue': eigenList.get(2),
    'VarianceExplained': varianceList.get(2)
  }),
  ee.Feature(null, {
    'Component': 'PC4',
    'ComponentNum': 4,
    'Eigenvalue': eigenList.get(3),
    'VarianceExplained': varianceList.get(3)
  })
]);


// Scree Plot - Autovalori
var screeChart1 = ui.Chart.feature.byFeature(screeFeatures, 'ComponentNum', ['Eigenvalue'])
  .setChartType('LineChart')
  .setOptions({
    title: 'Scree Plot - Autovalori',
    hAxis: {title: 'Componente Principale'},
    vAxis: {title: 'Autovalore'},
    lineWidth: 2,
    pointSize: 6,
    series: {
      0: {color: 'darkblue', pointShape: 'circle'}
    },
    legend: {position: 'none'}
  });
print(screeChart1);

// Scree Plot - Varianza spiegata (%)
var screeChart2 = ui.Chart.feature.byFeature(screeFeatures, 'ComponentNum', ['VarianceExplained'])
  .setChartType('LineChart')
  .setOptions({
    title: 'Scree Plot - Varianza Spiegata',
    hAxis: {title: 'Componente Principale'},
    vAxis: {title: 'Proporzione di Varianza Spiegata', format: 'percent'},
    lineWidth: 2,
    pointSize: 6,
    series: {
      0: {color: 'darkgreen', pointShape: 'circle'}
    },
    legend: {position: 'none'}
  });
print(screeChart2);

// Varianza cumulativa
var cumVar1 = varianceList.get(0);
var cumVar2 = ee.Number(cumVar1).add(varianceList.get(1));
var cumVar3 = ee.Number(cumVar2).add(varianceList.get(2));
var cumVar4 = ee.Number(cumVar3).add(varianceList.get(3));

var cumulativeFeatures = ee.FeatureCollection([
  ee.Feature(null, {'ComponentNum': 1, 'CumulativeVariance': cumVar1}),
  ee.Feature(null, {'ComponentNum': 2, 'CumulativeVariance': cumVar2}),
  ee.Feature(null, {'ComponentNum': 3, 'CumulativeVariance': cumVar3}),
  ee.Feature(null, {'ComponentNum': 4, 'CumulativeVariance': cumVar4})
]);

// Grafico varianza cumulativa
var cumulativeChart = ui.Chart.feature.byFeature(cumulativeFeatures, 'ComponentNum', ['CumulativeVariance'])
  .setChartType('LineChart')
  .setOptions({
    title: 'Varianza Cumulativa',
    hAxis: {title: 'Componente Principale'},
    vAxis: {title: 'Varianza Cumulativa', format: 'percent'},
    lineWidth: 2,
    pointSize: 6,
    series: {
      0: {color: 'purple', pointShape: 'circle'}
    },
    legend: {position: 'none'}
  });
print(cumulativeChart);
