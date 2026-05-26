#!/bin/bash
# pi-pi stacking validation
# Benchmark cases from Supplementary Table S2 (cases 1-11)
# Run from this directory: bash commands.sh

ncip 1eve.pdb -l E20 -c A -r 2001 --txt report_1eve.txt
ncip 1l8b.pdb -l MGP -c A -r 1000 --txt report_1l8b.txt
ncip 2w26.pdb -l RIV -c A -r 1001 --txt report_2w26.txt
ncip 3i1y.pdb -l 33N -c A -r 999  --txt report_3i1y.txt
ncip 4tvj.pdb -l 09L -c A -r 601  --txt report_4tvj.txt
ncip 4bfq.pdb -l 083 -c B -r 301  --txt report_4bfq.txt
ncip 1n5r.pdb -l PRM -c A -r 951  --txt report_1n5r.txt
ncip 18gs.pdb -l GDN -c A -r 210  --txt report_18gs.txt
ncip 1dmw.pdb -l HBI -c A -r 700  --txt report_1dmw.txt
ncip 1tll.pdb -l FAD -c A -r 1452 --txt report_1tll.txt
ncip 1f20.pdb -l FAD -c A -r 1501 --txt report_1f20.txt
