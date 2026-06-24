# Zhang/Yoon Reference Curve Digitization Notes

This file records the provenance and approximate digitization rules for `zhang_yoon_reference_curves.csv`.

## Target Times

The current CSV digitizes the benchmark image-comparison times used by the local precipitation runner:

```text
13 min
18 min
118 min
```

## Zhang 2010

Source:

```text
Zhang Table 1: experiment 3 = 25 mM with Omega_c/Omega_v = 4.6/3.9.
Zhang Supporting Information Figure S3(b): selected pore space at Omega_c/Omega_v = 4.6/3.9.
```

Regions:

```text
zhang_upgradient   -> Upgradient curve in SI Figure S3(b)
zhang_middle       -> middle curve in SI Figure S3(b)
zhang_downgradient -> Downgradient curve in SI Figure S3(b)
```

Calibration used for the local extracted image `zhang2010_si_page-4.png`:

```text
Panel S3(b) x-axis: 0 to 5 h
Panel S3(b) y-axis: 0 to 1.2E+04 pixels
Approximate plot bounds: x = 542 to 1186 px; y = 1411 to 937 px
```

Zhang reports area as pixels, not square centimeters. Therefore `precipitated_area_norm` is digitized as:

```text
digitized_pixel_value / 12000
```

and `precipitated_area_cm2` is left as `NaN`.

Approximate digitized pixel values:

```text
region                13 min   18 min   118 min
zhang_upgradient       4600     3850      2700
zhang_middle           6400     6900      2500
zhang_downgradient     8800    10500      4800
```

## Yoon 2012

Sources:

```text
Yoon Figure 4(a): Case 1 entire-domain precipitate area.
Yoon Figure 7(a): Case 5 entire-domain and first-pore precipitate area.
```

Yoon reports area in square micrometers. The CSV converts to square centimeters using:

```text
1 um2 = 1e-8 cm2
```

Approximate digitized values before unit conversion:

```text
source/region               13 min   18 min   118 min
case_1 entire_domain         36200    36400     36600 um2
case_5 entire_domain         34500    30500     17500 um2
case_5 first_pore            11000     9500      5600 um2
```

`precipitated_area_norm` is left as `NaN` for Yoon rows. `precip_CompareZhangYoonBenchmark` fills missing normalized values from `precipitated_area_cm2` per source/case/region using that group's maximum absolute area.

## Accuracy Boundary

These rows are visual digitizations from rasterized local figure images, not publisher tabular data. They are appropriate for a first benchmark overlay and source-bounded sanity check. They should not be presented as exact experimental data without either manual WebPlotDigitizer verification or source tabular data.
