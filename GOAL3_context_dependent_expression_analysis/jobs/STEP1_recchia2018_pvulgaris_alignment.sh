#!/bin/sh

#The purpose of this script is to match the reads from Recchia et al 2018 (https://doi.org/10.3389/fmicb.2018.01339) to the phaseolus vulgarum reference on NCBI. Judging from the read lengths (not stated in the paper). The reads are 100bp paired ends.

#The reads deposited in SRA are filtered for quality. The assumed working directory of this script is the jobs folder.

cd ../star_genomes
mkdir phaseolus_vulgaris_100bp

STAR --runMode genomeGenerate \
	--genomeDir phaseolus_vulgaris_100bp \
	--genomeFastaFiles ../raw_data/ncbi_genome_references/GCF_000499845.1_PhaVulg1_0_genomic.fna \
	--sjdbGTFfile ../raw_data/ncbi_genome_references/GCF_000499845.1_PhaVulg1_0_genomic.gtf \
	--sjdbOverhang 99 \
	--runThreadN 3 \
	--outFileNamePrefix phaseoulus_vulgaris_100b
