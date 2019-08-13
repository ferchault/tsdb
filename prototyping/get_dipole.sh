#!/bin/bash

#this script gives the dipole magnutude in Debye

logfile=$1

grep "Magnitude (Debye)      : " $logfile | tail -1 | awk '{print $4}' 
