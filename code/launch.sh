#!/bin/bash
# launch.sh
if [ -d "software/deps_r" ] 
then
    Rscript Installation_packages.R
    Rscript covid_network.R
else
    sudo bash setup.sh
fi