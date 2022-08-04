#!/bin/bash

# Wrapper for Phase I, II, III

# one gene at a time; usage: plink_wrapper.sh -c

while getopts "c" option
do 
	case $option in
		c) plot_overlay=true;; # if user wants make overlay plots; optional
	esac
done

gene_list="$(awk '{if(NR>1) print $2}' phenotype_for_chosen_genes | sort | uniq)"

for gene in $gene_list
do
	# For each gene:
	# 1. extract_phenotypes.sh -a $gene
	./extract_phenotypes.sh -a $gene

	# 2. run_plink.sh -a $gene
	./run_plink.sh -a $gene

	# 3. make csv files: csv_plots.sh
	#  make overlayed plots with -c
	if [ $plot_overlay ]; then
		./csv_plots.sh -a $gene -c # make plot input csv files and plots
	else
		./csv_plots.sh -a $gene # only make plot input csv files
	fi
done

echo "Finished!"

ls */*DISAGREED*
