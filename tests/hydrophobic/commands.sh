#!/bin/bash
# Hydrophobic contact validation
# Benchmark cases from Supplementary Table S2 (cases 1-7)
# Run from this directory: bash commands.sh

ncip 5c5p.pdb -l 0E0 -c A -r 1202 --txt report_5c5p.txt
ncip 2fra.pdb -l CRV -c A -r 999  --txt report_2fra.txt
ncip 1j4h.pdb -l SUB -c A -r 201  --txt report_1j4h.txt
ncip 2onz.pdb -l TMJ -c A -r 2001 --txt report_2onz.txt
ncip 3kwj.pdb -l 23Q -c A -r 1    --txt report_3kwj.txt
ncip 1m17.pdb -l AQ4 -c A -r 999  --txt report_1m17.txt
ncip 2yxj.pdb -l N3C -c A -r 1001 --txt report_2yxj.txt
