#!/bin/bash
# Water bridge validation
# Benchmark cases from Supplementary Table S2 (cases 1-7)
# Run from this directory: bash commands.sh

ncip 1hpx.pdb -l KNI -c B -r 900  --txt report_1hpx.txt
ncip 1z6f.pdb -l BO9 -c A -r 401  --txt report_1z6f.txt
ncip 5j1x.pdb -l AR5 -c A -r 501  --txt report_5j1x.txt
ncip 1b56.pdb -l PLM -c A -r 136  --txt report_1b56.txt
ncip 1ey3.pdb -l DAK -c A -r 500  --txt report_1ey3.txt
ncip 1kt4.pdb -l RTL -c A -r 184  --txt report_1kt4.txt
ncip 1tt2.pdb -l THP -c A -r 501  --txt report_1tt2.txt
