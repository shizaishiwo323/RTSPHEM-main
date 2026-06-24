# Zhang/Yoon Reference Curves

`zhang_yoon_reference_curves.csv` stores source-bounded approximate digitized reference rows. See `digitization_notes.md` for axis calibration, unit conversion, and accuracy limits.

Do not add placeholder or synthetic values to this file. Rows should be added only after digitizing Zhang 2010 experiment data and Yoon 2012 Case 1 / Case 5 curves from the source figures.

Required columns:

```text
source,case,region,time_min,precipitated_area_norm,precipitated_area_cm2,note
```

Allowed source/case/region combinations:

```text
Zhang2010 / experiment_25mM -> zhang_upgradient, zhang_middle, zhang_downgradient
Yoon2012  / case_1          -> entire_domain, first_pore, first_three_pores
Yoon2012  / case_5          -> entire_domain, first_pore, first_three_pores
```

Zhang 25 mM source basis:

```text
Zhang Table 1: experiment 3 = 25 mM, Omega_c/Omega_v = 4.6/3.9.
Zhang Supporting Information Figure S3(b): selected pore space at Omega_c/Omega_v = 4.6/3.9.
```

Yoon source basis:

```text
Yoon Figure 4(a): Case 1 and experiment curves for the entire numerical domain.
Yoon Figure 7(a): Case 5 and experiment curves for entire domain and first pore.
```

For Zhang pixel-based rows, record the digitized normalized value in `precipitated_area_norm` and leave `precipitated_area_cm2` as `NaN`. For Yoon square-micrometer rows, convert to square centimeters in `precipitated_area_cm2` and leave `precipitated_area_norm` as `NaN`; `precip_CompareZhangYoonBenchmark` derives normalized values downstream by source/case/region when needed.
