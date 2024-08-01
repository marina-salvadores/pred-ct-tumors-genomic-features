#!/usr/bin/bash

for filename in ../../results/ssm_remove/ssm_conOG/*; do
    echo $filename
    #Rscript RMD_calculation.R -i $filename -o '../../results/ssm_remove/RMD/RMD_'${filename:36:200} -f -c ../../data/raw_data/windows_ranges_coverage_alignability36.csv -t 100000 -r mean
    #Rscript MS96_calculation.R -i $filename -o '../../results/ssm_remove/MS96/MS96_'${filename:36:200}
    Rscript OG_calculation.R -i $filename -o '../../results/ssm_remove/OG/OG_'${filename:36:200} -p -s

done
