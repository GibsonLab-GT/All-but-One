# AbO-Conditional-Analysis
A pipeline to perform "all-but-one" conditional analysis with the objective of isolating eQTL peaks by reducing linkage disequilibrium interference.


----------
OBJECTIVE:
----------
The objective of this pipeline is to perform "all-but-one" conditional analysis by isolating sequential eQTL peaks and reducing effects of linkage disequilibirum (LD).

-------------------
STEPS & HOW TO RUN:
-------------------

I.) Input:

	Genotype data for each gene of interest in the form of PLINK files located in a "genotype" directory. Any phenotype data must be called "phenotype_for_chosen"genes" and structured as:
	
	chr
	gene
	probe/misc
	quantitative phenotypic value
	
	These scripts are required to run the pipeline, along with PLINK. PLINKv1.9 can be accessed and downloaded here: https://www.cog-genomics.org/plink/
	
	i.) plink_wrapper.sh
	ii.) extract_phenotypes.sh
	iii.) transpose.sh
	iv.) run_plink.sh
	v.) csv_plots.sh
	vi.) LZ_overlay_plots.R --> requires R package libraries "readr" and "stringi"
	vii.) csv_plots_No_Cond_Analysis.sh

II.) Automated Steps:

	For each gene listed in the "phenotype_for_chosen_genes" file:
		Generate a new directory
		copy the genotype data into this directory
		
		For each probe associated with this gene
			create a PLINK compatible phenotype file by subsetting the "phenotype_for_chosen_genes" file
			perform stepwise linear association analysis with PLINK, noting each eSNP meeting significance threshold
			
			For each eQTL peak discovered
				use the forward stepwise "eSNP list" to perform "All but One" conditional analysis, noting any potential eSNP disagreements
				compute LD R^2 between the eSNP versus all other SNPs in this region
				generate custom locuszoom-like plots

III.) Running the script:

	i.) version 1: run automated pipeline to all genes in "phenotype_for_chosen_genes" file:

		./plink_wrapper.sh -c
	
	ii.) version 2: run pipeline for ONE gene in "phenotype_for_chosen_genes" file (all probes for this gene run automatically):
		
		Example for CISD1:
			./extract_phenotypes.sh -a CISD1
			./run_plink.sh -a CISD1
			./csv_plots.sh -a CISD1 -c  

IV.) Plotting Forward Stepwise locuszoom-like plots (NOT AbO conditional analysis locuszoom-like plots)

	Example for CISD1:
		./csv_plots_No_Cond_Analysis.sh -a CISD1
