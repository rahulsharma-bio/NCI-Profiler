#!/bin/bash
# Metal coordination validation
# Benchmark cases from Supplementary Table S2 (cases 1-13)
# Run from this directory: bash commands.sh

ncip 4hhb.pdb -l HEM -c A -r 142  --txt report_4hhb.txt
ncip 2bb6.pdb -l B12 -c A -r 0   --txt report_2bb6.txt
ncip 1a42.pdb -l BZU -c A -r 555  --txt report_1a42.txt
ncip 1azm.pdb -l AZM -c A -r 262  --txt report_1azm.txt
ncip 3hs4.pdb -l AZM -c A -r 701  --txt report_3hs4.txt
ncip 1jk3.pdb -l BAT -c A -r 900  --txt report_1jk3.txt
ncip 1o86.pdb -l LPR -c A -r 702  --txt report_1o86.txt
ncip 1t64.pdb -l TSN -c A -r 386  --txt report_1t64.txt
ncip 1w22.pdb -l NHB -c A -r 1378 --txt report_1w22.txt
ncip 3oya.pdb -l RLT -c A -r 398  --txt report_3oya.txt
ncip 3thh.pdb -l ABH -c A -r 700  --txt report_3thh.txt
ncip 5yxc.pdb -l CIT -c A -r 302  --txt report_5yxc.txt
ncip 6dzq.pdb -l HJP -c A -r 203  --txt report_6dzq.txt
