#!/bin/bash
# Cation-π interaction validation
# Benchmark cases from Supplementary Table S2 (cases 1-9)
# Run from this directory: bash commands.sh

ncip 2jkh.pdb -l BI7 -c A -r 1244 --txt report_2jkh.txt
ncip 2bok.pdb -l 784 -c A -r 1244 --txt report_2bok.txt
ncip 1acj.pdb -l THA -c A -r 999  --txt report_1acj.txt
ncip 2y5h.pdb -l Y5H -c A -r 1244 --txt report_2y5h.txt
ncip 1uv6.pdb -l CCE -c C -r 1206 --txt report_1uv6.txt
ncip 1uw6.pdb -l NCT -c A -r 1208 --txt report_1uw6.txt
ncip 4mqs.pdb -l IXO -c A -r 501  --txt report_4mqs.txt
ncip 3gwv.pdb -l RFX -c A -r 801  --txt report_3gwv.txt
ncip 2b4l.pdb -l BET -c A -r 273  --txt report_2b4l.txt
