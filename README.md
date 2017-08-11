# GxE_scripts

The *stability_associations* directory contains scripts used to perform Finlay-Wilkinson regression on G2F 2014 phenotypic data, run GWAS on the stability parameters using GAPIT, and characterize the strongest association results.

**1_pheno_prep.R**: Compiles check mean data, aggregates hybrid phenotypes across locations, and regresses each hybrid on the check means.  Makes lots of plots.

**2_merge_phenos_forGWAS.R**: Cleans regression parameters, plots cleaned data, writes out clean data to file.

**3_format_genos_forGAPIT.R**: Formats hybrid genotypes into format compatible with GAPIT.

**4_run_GAPIT.R**: runs GAPIT with cleaned regression parameters as the response variables.

**5_calc_LD_plink.sh**: calculates LD between markers within a 5kb window of each other.

**6_analyze_GAPIT_results.R**: creates manhattan and QQ plots of GAPIT results and writes the most significant SNPs out to a new file.

**7_findClosest_Genes_withDirection.R**: finds the distance to the closest gene for each SNP from the previous script. Performs and exact binomial test against the null distribution.  Makes a plot.

The *nucleotide_diversity* directory contains scripts to perform nucleotide diversity analysis within 30 temperate and 30 tropical inbred lines.

**1_get_BAM_files.sh**: Pulls the sequence files for the inbreds of interest down from Cyverse

**2_create_merge_sample_files.sh**: Some inbreds were sampled and sequenced more than once - this lists the sample files for each inbred line.

**3_merge_samples.sh**: Merges the sequence files for inbreds with >1 sequence files

**4_get_high_and_low_Fst_windows.R**: Gets the bounds of high and low Fst windows

**5_get_SAFs.sh**: Calculates the site allele frequencies within the temperate and tropical groups

**6_makeSFS.sh**: Uses the SAF files to create a site frequency spectrum for temperate and tropical groups

**7_nucleotide_diversity.sh**: Calculates nucleotide diversity within the temperate and tropical groups

**8_merge_nd.sh**: Pulls the nucleotide diversity results together into a usable format

**9_window_nd.R**: Nucleotide diversity statistics are on a per-site basis -- this calculates per-window averages

**10_merge_windows.sh**: Pulls all the per-window nucleotide diversity stats together into one file

**11_plot_nd.R**: Creates the figure used in the publication
