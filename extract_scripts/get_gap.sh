#!/bin/bash

#this script gives the homo lumo and gap energies in Hartree

logfile=$1

ne=$(grep " Number of Electrons    NEL             ...." $logfile | head -1 | awk '{print $6}')

lumo=$(( $ne / 2 | bc -l))

homo=$(( $lumo -1 | bc -l))

E_lumo=$(grep " $lumo   0.0000   " $logfile | tail -1 | awk '{print $3}')

E_homo=$(grep " $homo   2.0000   " $logfile | tail -1 | awk '{print $3}')

gap=$( echo $E_lumo- $E_homo | bc -l)

echo "$E_homo  $E_lumo  $gap"
