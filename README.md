# GENOMIC FEATURES MATRIX CALCULATION

## 1. Regional mutation density (RMD)
__RMD_calculation.R__ - script that calculates a RMD matrix from a single somatic mutations file.

Run script:
```bash
Rscript RMD_calculation.R -i name_input_ssm_file -o name_output -f -c bp_cov_per_window_file -t 100000 -r mean
```

Parameters needed:
- __"-i, --input"__ ssm file (csv containing chr, start, end, sample_id, cancer_type)
- __"-o, --output"__ output path and name for the RMD matrix
- __"-f, --norm_by_factor"__ add -f when you want to normalize by factor. For this 2 additional parameters are required:
  - __"-c, --cov_file"__ file containing the windows coverage for the ssm data type (alignability of simuWXS)
  - __"-t, --cov_threshold"__ threshold for removing windows lower than the specify bp coverage
- __"-r, --norm_by_row"__ normalize by row. Options: 'mean' or 'sum'



## 2. Mutation spectrum (MS96)
__MS96_calculation.R__ - script that calculates a MS96 matrix from a single somatic mutations file.

Run script:
```bash
Rscript MS96_calculation.R -i name_input_ssm_file -o name_output
```

Parameters needed:
- __"-i, --input"__ ssm file (csv containing chr, start, end, sample_id, cancer_type)
- __"-o, --output"__ output path and name for the MS96 matrix



## 3. Presence/absence of mutations in OG and TSG (OG)
__OG_calculation.R__ - script that calculates a OG matrix from a single somatic mutations file.

Run script:
```bash
Rscript OG_calculation.R -i name_input_ssm_file -o name_output
```

Parameters needed:
- __"-i, --input"__ ssm file (csv containing chr, start, end, sample_id, cancer_type)
- __"-o, --output"__ output path and name for the MS96 matrix
- __"-p, --pathways"__ add -p to calculate pathways features and add them to og matrix
- __"-s, --hotspots"__ add -h to calculate hotspots features and add them to og matrix


## 4. Total number of mutations (SUM)
__SUM_calculation.R__ - script that calculates the SUM matrix from a single somatic mutations file.

Run script:
```bash
Rscript SUM_calculation.R -i name_input_ssm_file -o name_output
```

Parameters needed:
- __"-i, --input"__ ssm file (csv containing chr, start, end, sample_id, cancer_type)
- __"-o, --output"__ output path and name for the SUM matrix


## 5. Copy number alterations (CNAs)
__CNA_calculation.R__ - script that calculates the CNAs matrix from a folder of CNAs files (per cantype).

Run script:
```bash
Rscript CNAs_calculation.R -i name_input_file -o name_output
```

Parameters needed:
- __"-i, --input"__ any file that contains the sample_id and cancer_type of the samples we are interested in
- __"-o, --output"__ output path and name for the CNAs matrix
