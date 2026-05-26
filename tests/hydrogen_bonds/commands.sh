#!/bin/bash
# Hydrogen bond validation
# Benchmark cases from Supplementary Table S2 (cases 1-11)
# Run from this directory: bash commands.sh

ncip 2zc9.pdb -l 22U -c H -r 1501 --txt report_2zc9.txt
ncip 1pxp.pdb -l CK8 -c A -r 500  --txt report_1pxp.txt
ncip 3ehx.pdb -l BDL -c A -r 0    --txt report_3ehx.txt
ncip 2c3l.pdb -l IDZ -c A -r 1274 --txt report_2c3l.txt
ncip 2c3k.pdb -l ABO -c A -r 1271 --txt report_2c3k.txt
ncip 1owe.pdb -l 675 -c A -r 1001 --txt report_1owe.txt
ncip 1hvr.pdb -l XK2 -c A -r 263  --txt report_1hvr.txt
ncip 1fkb.pdb -l RAP -c A -r 108  --txt report_1fkb.txt
ncip 3ert.pdb -l OHT -c A -r 600  --txt report_3ert.txt
ncip 4hjo.pdb -l AQ4 -c A -r 1001 --txt report_4hjo.txt
ncip 4ekl.pdb -l 0RF -c A -r 501  --txt report_4ekl.txt
