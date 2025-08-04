#!/bin/sh

#The goal in this script is to perform MCSanX_h on all the species used in our OrthoFinder analysis.
#Downstream, we will remove taxa that had too many scaffolds because synteny analysis is very difficult when chromosomes are too split up.
#The assumed working directory is the jobs folder

#get path to MCScanX_h
mcscan=/home/damian/nonstandard_programs/MCScanx

#make list of species that were included in OrthoFinder analysis
specs="AM_aegilops_tauschii AM_ananas_comosus AM_asparagus_officinalis AM_brachypodium_distachyon AM_cajanus_cajan AM_capsicum_annuum AM_carica_papaya AM_cicer_arietinum AM_cucumis_melo AM_daucus_carota AM_fragaria_vesca AM_gossypium_arboreum AM_helianthus_annuus AM_hevea_brasiliensis AM_lactuca_sativa AM_manihot_esculenta AM_medicago_truncatula AM_musa_acuminata AM_nicotiana_tabacum AM_olea_europaea AM_phaseolus_vulgaris AM_phoenix_dactylifera AM_prunus_persica AM_ricinus_communis AM_sesamum_indicum AM_setaria_italica AM_solanum_lycopersicum AM_sorghum_bicolor AM_theobroma_cacao AM_vigna_radiata AM_vitis_vinifera AM_zea_mays"

#move into processed_data directory for for loop to be able to go back and forth between species
cd ../processed_data

for spec in $specs
	do
		cd $spec
		$mcscan/MCScanX_h ${spec}_mcscan_input -b 1
		cd ..
	done

cd ../jobs
