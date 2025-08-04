#!/bin/sh

#The purpose of this script is to download the SRA data associated with Recchia et al 2018.
#The assumed working directory is the jobs folder.

cd ../raw_data
mkdir recchia2018_sra
cd recchia2018_sra

prefetch --option-file ../../metadata/recchia2018_SRR_Acc_List.txt
