#!/bin/bash
# XWH (halogen-water-hydrogen) bridge validation
# Benchmark cases from Supplementary Table S2 (cases 1-10)
# Run from this directory: bash commands.sh

ncip 2uy3.pdb -l H33 -c A -r 1311 --txt report_2uy3.txt
ncip 2viv.pdb -l VG2 -c A -r 1247 --txt report_2viv.txt
ncip 1g3m.pdb -l PCQ -c A -r 712  --txt report_1g3m.txt
ncip 1g5f.pdb -l DCE -c A -r 601  --txt report_1g5f.txt
ncip 1uv5.pdb -l BRW -c A -r 1383 --txt report_1uv5.txt
ncip 1n8u.pdb -l BDD -c A -r 411  --txt report_1n8u.txt
ncip 1sn5.pdb -l T3  -c C -r 601  --txt report_1sn5.txt
ncip 3e3d.pdb -l I3C -c A -r 130  --txt report_3e3d.txt
ncip 2b1p.pdb -l AIZ -c A -r 501  --txt report_2b1p.txt
ncip 1wcc.pdb -l CIG -c A -r 1299 --txt report_1wcc.txt
