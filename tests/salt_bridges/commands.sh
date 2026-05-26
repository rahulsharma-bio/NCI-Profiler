#!/bin/bash
# Salt bridge validation
# Benchmark cases from Supplementary Table S2 (cases 1-9)
# Note: 7q2n.pdb and 7q19.pdb were converted from mmCIF using gemmi
# (original PDB format not available from RCSB)
# Run from this directory: bash commands.sh

ncip 4jrv.pdb -l KJV -c A -r 1001 --txt report_4jrv.txt
ncip 2z7k.pdb -l BGU -c A -r 400  --txt report_2z7k.txt
ncip 1bji.pdb -l DPC -c A -r 479  --txt report_1bji.txt
ncip 1ppc.pdb -l MID -c E -r 5    --txt report_1ppc.txt
ncip 3pbl.pdb -l ETQ -c A -r 1200 --txt report_3pbl.txt
ncip 7q2n.pdb -l DSM -c A -r 201  --txt report_7q2n.txt
ncip 7q19.pdb -l DSM -c A -r 201  --txt report_7q19.txt
ncip 5mu6.pdb -l KFK -c A -r 505  --txt report_5mu6.txt
ncip 3odu.pdb -l ITD -c A -r 1500 --txt report_3odu.txt
