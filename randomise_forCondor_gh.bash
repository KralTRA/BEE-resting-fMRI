#!/bin/sh

inputFile=$1
design=${inputFile}

datPath="/path/"

randomise -i ${datPath}${inputFile}_4D.nii.gz -o ${datPath}${inputFile} -d ${datPath}${design}.mat -t ${datPath}${design}.con -n 5000 -T