#! /bin/bash

### Stephen JD Kay, University of Regina
### 25/05/21
### stephen.kay@uregina.ca
### A script to process a bunch of .lund DEMPgen outputs through EIC smear.
### For now, paths are hardcoded
### /data is simply the bound directory containing the data files on my machine

i=1
while [[ $i -le 100 ]]; do

    root -l<<EOF
gSystem->Load("libeicsmear.so");
BuildTree("/data/DEMP/LundFiles/5_100_eic_smear_feed_in_21_05_21/1B_Event_Files/eic_DEMPGen_5on100_1000000000_${i}.lund", ".", -1, "log_${i}.txt")
.q
EOF
    i=$(( $i + 1 ))

done
