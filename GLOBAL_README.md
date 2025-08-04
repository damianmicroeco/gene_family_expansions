## Preface
The goal of this project is to determine whether gene family expansions are integral to context-dependent regulation of arbuscular mycorrhizal symbiosis. All files to replicate our analyese, with the exception of the HapMap SNP data, are available in this folder. For the computational environments, please find these at the end of this README file.

Please not that for this github repository, I am only including the scripts because a lot of the data is too big for GitHub's size requirements.
However, you can find this data on Zenodo. Please find the citation below:
Hernandez, D. (2024). Gene family expansions underpin context-dependency of the oldest mycorrhizal symbiosis [Data set]. Zenodo. https://doi.org/10.5281/zenodo.14340319

## General Structure
This project folder is separated into 5 "goals" folders and 1 plotting folder. All the "goals" folders are self-contained and will not reference anything outside itself. The plotting script is the only folder that will occassionally search for files outside of its own folder. All scripts assume the working directory is the jobs folder within each "goals"/plotting folder. 

The general structure of most folders are as follows:
* raw data: the original data used for analyses without any edits. So, these are data that were directly downloaded from their source (e.g., SRA, manuscripts, ENSEMBL Plants, etc).
* processed_data: data that have been modified by the scripts.
* metadata: data that may be original or edited by a computational script. These are data that mostly provide information about identity.
* jobs: the folder from which scripts are run. The order in which scripts are run are stated in the name of the file. For example, the script prefaced with "STEP1" is run first and the one prefaced with "STEP2" is run second.

Certain folders may have some other folders needed for a specific or slightly modified structure due to the algorithm that was used. These are outlined below in the more detailed explanations about the different folders.

## GOAL1_isoform_selection
The main goal in this folder is to filter the RefSeq protein data for the longest isoform of each gene. This step is necessary for input into OrthoFinder because, without selecting one isoform for each gene, then orthgroup sizes are inflated by the number of isoforms in that orthogroup rather than the number of genes. The assumed working directory is the root of this GOAL folder.

Folders and Files:
* input_data: data for input to scripts. This is equivalent to the "raw_data" folders from the general structure.
    * genome_metadata: tab-delimited files of NCBI RefSeq information for all plant species mapping RefSeq protein IDs to GeneIDs at the time of download.
    * refseq_protein: FASTA protein files from RefSeq for all plant species in analysis
* misc: log files of proteins that could not be matched between the genome metadata and RefSeq FASTA files.
    * missing_proteins: the lists of missing proteins for each species
    * missing_proteins_metadata.txt: log file summarizing how many proteins were found in each species.
* processed_data
    * longest_protein: Lists of RefSeq IDs that were identified as the longest isoforms for their GeneID
    * refseq_longest_select: FASTA files of the selected isoforms in the longest_protein folder.
    * refseq_longest_select: Manually selected FASTA files from refseq_longest_select for species that match selection criteria
* scripts
    * python_scripts
        * STEP1_max_isoform_selection.py: script identifying the longest isoforms of each GeneID
        * STEP2_max_isoform_FAA_create.py: script generating the FASTA files of the longest isoforms of each GeneID for input into OrthoFinder.

## GOAL2_orthogroup_analysis
The main goal of this folder is to store the output from the commandline OrthoFinder function which was provided with default parameters and the refseq_longest folder in GOAL1 as input data. This folder only contains the output of the OrthoFinder analysis rather than running it because it is a one-line command (e.g. `orthofinder -f refseq_longest_select/`)

Folders and Files:
* Orthologues_Sep12: Intermediate files produced by OrthoFinder to build the orthogroup classifications
    * Duplications.csv
    * OrthologuesStats_one-to-one.csv
    * Gene_Trees
    * OrthologuesStats_Totals.csv
    * Orthologues
    * Recon_Gene_Trees
    * OrthologuesStats_many-to-many.csv
    * SpeciesTree_rooted.txt
    * OrthologuesStats_many-to-one.csv
    * SpeciesTree_rooted_node_labels.txt
    * OrthologuesStats_one-to-many.csv
    * WorkingDirectory
* Orthogroups.GeneCount.csv: tab-separated text file containing number of genes in each orthogroup in each species. rows are orthogroups, columns are species.
* Orthogroups_UnassignedGenes.csv: tab-separated text files with orthogroups as rows and species as columns. Cells are the genes in that species that belong to that orthogroup, but these are genes that could not be assigned to a real orthogroup.
* Orthogroups.tsv: tab-separated text files with orthogroups as rows and species as columns. Cells are the genes in that species that belong to that orthogroup.
* Orthogroups.txt: legacy file format of the OrthoFinder program repeating the data in Orthogroups.tsv, but in the OrthoMCL output format.
* SingleCopyOrthogroups.txt: orthogroups that have only one gene in each species.
* orthogroups_mannwhitney_sig.csv: Results of Mann-Whitney U comparisons between orthogroup sizes in AM and NM plants.
* Statistics_Overall.csv: tab-separated log-file of orthogroup classifications for all species together
* Orthogroups_SpeciesOverlaps.csv: tab-separated text file containing matrix of the number of orthogroups shared by each species pair.
* Statistics_PerSpecies.csv: tab-separated log-file of orthogroup classifications per species

## GOAL3_context_dependent_expression_analysis
The main goal in this folder is to conduct context-dependent expression analyses of genes in factorial experiments of plants grown with/without mycorrhizal fungi and with/without an environmental stressor. Some preparatory work was sometimes needed to, for example, align raw data from SRA to RefSeq genes, match ENSEMBL IDs to RefSeq IDs, etc.

Folders and Files:
* 20200324_genome_analyses: This folder contains data from GOAL1 and GOAL2 just to keep this folder contained.
* blast_db: folder to store constructed BLAST databases for reference alignments
* jobs: folder containing analyses. These are run in the orders they are prefaced. There are four separate analyses so which STEP2 to run first does not matter as long as all the STEP1 jobs have been run.
    * STEP0-5: this is to download the raw sequencing data from Recchia et al 2018.
    * STEP1: These are to match RefSeq IDs to ENSEMBL IDs or to align raw sequencing data through STAR aligner to a reference genome.
    * STEP2: Context-dependent expression analyses that compares DESeq2 model containing interaction terms with those that do not and then assesses if the orthogroups with more genes have more context-dependent expression. The observed values in the AM-expanded orthogroups are compared to 10,000 random subsamples across the whole genome to determine signifigance.
* metadata: metadata files mostly containing identity information for genes
    * calabrese2019_metadata.csv: file containing information about which samples belonged to which treatment group in Calabrese et al 2019
    * garcia2017_metadata.csv: file containing information about which samples belonged to which treatment group in Garcia et al 2017
    * manihot_esculenta_proteins_441_277153.csv: RefSeq table for Manihot esculenta
    * Medicago_truncatula.MedtrA17_4.0.51.ena.tsv: ENSEMBL Plants table for Medicago truncatula.
    * orthogroup_identified_AM_medicago_truncatula_genome_metadata.tsv: RefSeq table for Medicago truncatula with orthogroup identity
    * orthogroup_mannwhitney.csv: Results of Mann-Whitney U comparisons between orthogroup sizes in AM and NM plants.
    * Phaseolus_vulgaris.PhaVulg1_0.51.entrez.tsv: ENTREZ table for Phaseolus vulgaris
    * Populus_trichocarpa.Pop_tri_v3.51.ena.tsv: ENSEMBL Plants table for Populus trichocarpa.
    * recchia2018_metadata.csv: file containing information about which samples belonged to which treatment group in Recchia et al 2018. Directly downloaded from NCBI
    * recchia2018_metadata_edited.csv: same as recchia2018_metadata.csv, but with an extra column to provide treatment information that is compatible with R.
    * recchia2018_SRR_Acc_List.txt: List of SRA accessions to download SRA data for Recchia et al 2018
    * ren2019_metadata.csv: file containing information about which samples belonged to which treatment group in Ren et al 2019
* processed_data
    * 20210706_medicago_blast_results.tsv: Results of BLAST matching between RefSeq and ENSEMBL for medicago truncatula
    * 20210706_medicago_blast_with_orthogroup.tsv: Results of BLAST matching between RefSeq and ENSEMBL for medicago truncatula and orthogroup identity
    * 20210713_sesbania_blast_results.tsv: Results of BLAST matching between Sesbania cannabina and Medicago truncatula
    * 20210713_sesbania_cannabina_blast_with_orthogroup.tsv:  Results of BLAST matching between Sesbania cannabina and Medicago truncatula with orthogroup identity
    * 20210718_populus_blast_with_orthogroup.tsv: Results of BLAST matching between Populus trichocarpa and Manihot esculenta with orthogroup identity
    * 20210823_phaseolus_blast_results.tsv: Results of BLAST matching between RefSeq and ENSEMBL for Phaseoulus vulgaris
    * 20210823_phaseolus_blast_with_orthogroup.tsv: Results of BLAST matching between RefSeq and ENSEMBL for Phaseolus vulgaris and orthogroup identity
    * am_orthogroup_contextdependentexpression: Folder containing DESeq2 results for all genes in each experiment
    * populus_blast_results.tsv: Results of BLAST matching between Populus trichocarpa and Manihot esculenta
    * Q-medicagoENSEMBL_S-medicagoREFSEQ: BLAST comparisons between ENSEMBL and RefSeq for Medicago truncatula.
    * Q-phaseolusENSEMBL_S-phaseoulusREFSEQ: BLAST comparisons between ENSEMBL and RefSeq for Phaseolus vulgaris
    * Q-populus_S-manihot: BLAST comparisons between Populus trichocarpa and Manihot esculenta
    * Q-ren2019_S-medicagoREFSEQ: BLAST comparisons between Sesbania cannabina and Medicago truncatula
    * recchia2018_pvulgaris_bam: BAM files for Phaseoulus vulgaris alignment for STAR aligner 
    * transcriptome_relationships: Files containing observed relationships between context-dependent expression in AM-expanded orthogroups and the relationships in the 10000 randomized subsamples.
* raw_data
    * calabrese_2019: raw data from Calabrese et al 2019
    * garcia2017_featurecounts_refseq.txt: raw data from Garcia et al 2017
    * Medicago_truncatula.MedtrA17_4.0.pep.all.fa: Medicago truncatula ENSEMBL data
    * medicago_truncatula_embl_seq: Medicago truncatula ENSEMBL data
    * ncbi_genome_references: genomes from NCBI
    * palakurty_2018: raw data from Palakurty et al 2018
    * palakurty_2018_rawcounts.csv: raw data from Palakurty et al 2018
    * Phaseolus_vulgaris.PhaVulg1_0.pep.all.fa: Phaseolus vulgaris ENSEMBL data
    * phaseolus_vulgaris_embl_seq: Phaseolus vulgaris ENSEMBL data
    * Populus_trichocarpa.Pop_tri_v3.pep.all.fa: Populus trichocarpa ENSEMBL data
    * populus_trichocarpa_embl_seq: Populus trichocarpa ENSEMBL data
    * recchia2018_sra: raw SRA data for recchia et al 2018
    * ren_2019: raw data for Ren et al 2019
    * ren2019_denovo_seq: Raw data for Ren et al 2019
* star_genomes: Files for STAR aligner for Phaseolus vulgaris

## GOAL4_GWAS
The main goal in this folder is to determine if larger AM-expanded gene families have more mycorrhizal-associated genetic variation with pod counts (i.e., fitness). We use the Medicago HapMap data. We do not include that data here because it is publicly available at https://medicago.legumeinfo.org/. If the organization no longer provides this data, please contact us and we will provide it. We also do not include the bcf to vcf conversions because they are exceedingly large and not highly modified from the original HapMap data. However, all the scripts to generate the processed files from the original HapMap data are in the `jobs` folder.

Folders and Files:
* jobs
    * ortho_size_snp_variation_pod.RData: RData used in size variation analyses from STEP4 onward (and in plotting script) so that analyses don't have to be re-run from the beginning for analysis. It is very time consuming.
    * STEP1_vcffiltering_gwas.sh: script to conver bcf file to vcf, filter for samples with fitness data, and to remove SNPs that do not meet inclusion criteria
    * STEP2_GWAS_Afkhami2021_preparatory.Rmd: script to prepare data for LFMMs in GWAS
    * STEP3_GWAS_pod.Rmd: script to conduct LFMMs for GWAS analysis
    * STEP4_GWAS_pod_SNP_mapping.Rmd: script to map SNP/gene/orthogroup identities to GWAS results
    * STEP5_ortho_size_snp_variation_pod.Rmd: script comparing observed relationship of percent sequence with SNP and orthogroup size in AM-expanded orthogroups to 10,000 random subsamples of the Medicago truncatula genome.
* metadata
    * mt4-0_jcvi: medicago genomic data with JCVI classifications
    * mt4-0_refseq: medicago genomic data with RefSeq classifications
    * 20210706_medicago_blast_with_orthogroup.tsv: BLAST results matching JCVI classifications to RefSeq with Orthogroup grouping
    * Medicago_truncatula.MedtrA17_4.0.52.refseq.tsv: Medicago truncatula genomic coordinates with RefSeq classifications
    * 20220407_vcffiltering_gwas_output_samplesnotinSNP.out: Log file of vcf filtering for samples not in HapMap SNP data.
    * Mt4.0_HapMap_README.pdf: README file from Medicago HapMap describing the data structure
    * Afkhami2021_data_by_genotype_dryad.csv: Raw performance data from Afkhami et al 2021. Downloaded from Dryad.
    * afkhami2021_investment_gwas.csv: Root-to-Shoot ratio data from Afkhami et al 2021 edited for input into R to match how GWAS analysis is structured
    * afkhami2021_pod_gwas.csv: Pod count data from Afkhami et al 2021 edited for input into R to match how GWAS analysis is structured
    * ncbi_dataset.zip: Raw data from NCBI for Medicago truncatula
    * hapmap_metadata.txt: Raw metadata from Medicago HapMap. Tab-delimited.
    * orthogroup_identified_AM_medicago_truncatula_genome_metadata.tsv: Genomic coordinates of Medicago truncatula genes with Orthogroup classification.
    * hapmap_metadata.xlsx: Raw metadata from Medicago HapMap in Excel format
    * orthogroup_mannwhitney.csv: Results of Mann-Whitney U comparisons between orthogroup sizes in AM and NM plants.
* processed_data
    * 20220407_GWAS_pod_output: LEA project output from GWAS for changes in pod investment
    * 20220407_GWAS_root2shoot_output: LEA project output from GWAS for changes in Root-to-Shoot Ratio
* raw_data
    * This folder is currently empty, but would normally house the SNP data directly from the Medicago Hapmap

## GOAL5_duplication_origin_analysis
The main goal in this folder is to determine if certain gene duplication events are enriched in AM-expanded orthogroups. Because the default gene duplication origin analysis function in MCScanX does not use orthogroups like we do, we recreated as much of the analysis as possible in R and leveraged the MCScanX_h function to do the whole-genome/segmental duplication analyses with orthogroup classifications.

Folders and Files:
* edited_images
    * 20231017_wgd_pvalue_summary: Summary of p-value results of duplication origin analysis.
    * 20231018_all_random_distributions_duplication_types: Violin plots of random distributions from duplication origin analysis
    * 20231019_wgd_pvalue_summary: Summary of p-value results of diplication origin analysis.
* jobs
    * duplicate_gene_classifier_all.RData: RData package holding all the results from the STEP2 function.
    * STEP2_duplicate_gene_classifier_all.Rmd: R recreation of MCScanX logic to quantify the duplication origins and compare them to random subsamples of the whole plant genome for all 24 species with high quality scaffolds.
    * STEP1_mcscanx_input_prep.ipynb: IPython notebook to prepare genomic data for input into MCScanX and/or the R recreation.
    * STEP2-5_MCScanX_h_iteration.sh: Shell script run in the middle of STEP2 to do the MCScanX_h function. Depending on your system, you can do this within the R script, but there are issues on Windows to access the Windows Subsystem for Linux as we do. So, we separated the MCScanX_h analysis into it's own shell script.
* metadata
    * genome_metadata: metadata containing genomic coordinates of RefSeq classifications for each protein accession
    * longest_protein: longest RefSeq isoforms for each gene
    * orthogroup_mannwhitney.csv: Results of Mann-Whitney U comparisons between orthogroup sizes in AM and NM plants.
    * Orthogroups.csv: tab-separated text files with orthogroups as rows and species as columns. Cells are the genes in that species that belong to that orthogroup.
    * Orthogroups.GeneCount.csv: tab-separated text file containing number of genes in each orthogroup in each species. rows are orthogroups, columns are species.
* processed_data: folders containing the distributions of duplication events, the random distributions, and all the input processing for use in MCScanX
* raw_data: placeholder folder. not used because all the data is in the metadata folder.
* temp_images: violin plots of duplications in random subsamples.

## PLOTTING
The main goal in this folder is to create the plots for the manuscript and some simple statistics.

Folders and Files:
* base plots
    * 20220914_fig3_transcriptome_edited.svg: base plot of context-dependent expression analysis.
    * 20220914_fig5_gwas_edited.svg: base plot of SNP variation analysis
    * 20231205_AMvNM_Heatmap_noplastid.svg: Heatmap distribution of orthogroups sizes across all plants
    * 20231130_duplication_types_plot.svg: base plot of duplication origin analysis
* customizations: Figures with touch-ups (e.g., editing axis labels, adding subfigure labels, making legends, etc) in InkScape
    * 20220623_SuppFig1_AMOrtho_Validation.png
    * 20231208_conceptualfigure.png
    * 20231211_Figure2_merged.png
    * 20231019_suppfig2_GWAS_PCA.png
    * 20231208_duplication_types_plot_edited.png
* intermediate_data
    * 20220517_heatmap_clustering_noplastid.csv: results of heatmap clustering in Figure 2
* jobs
    * 20231122_prettyplotting.Rmd: plotting script with simple statistics.