#!/bin/bash

# Wrapper for Phase I, II, III

# one gene at a time; usage: plink_wrapper.sh -c

while getopts "c" option
do 
	case $option in
		c) plot_overlay=true;; # for user to make plots; optional
	esac
done

gene_list="$(awk '{if(NR>1)print $2}' phenotype_for_chosen_genes | sort | uniq)"

for gene in $gene_list
do
	# For each gene and probe
	# Make csv files and plots
	elif [ $plot_overlay ]; then
		./csv_plots_No_Cond_Analysis.sh -a $gene -c # make plot input csv files and plots
	else
		./csv_plots_No_Cond_Analysis.sh -a $gene # only make plot input csv files, no plotting
	fi
done

echo "Finished!"
