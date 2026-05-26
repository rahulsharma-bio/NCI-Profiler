#!/bin/bash
# Halogen bond validation
# Benchmark cases from Supplementary Table S2 (cases 1-8)
# Run from this directory: bash commands.sh

ncip 5u00.pdb -l 7OV -c A -r 1001 --txt report_5u00.txt
ncip 1p5e.pdb -l TBS -c A -r 301  --txt report_1p5e.txt
ncip 3shz.pdb -l 5CO -c A -r 1    --txt report_3shz.txt
ncip 3sie.pdb -l 5BO -c A -r 1    --txt report_3sie.txt
ncip 3dn3.pdb -l IBF -c A -r 900  --txt report_3dn3.txt
ncip 4lbs.pdb -l 4O8 -c A -r 402  --txt report_4lbs.txt
ncip 4oew.pdb -l 5IO -c A -r 903  --txt report_4oew.txt
ncip 4x21.pdb -l 3WH -c A -r 501  --txt report_4x21.txt
