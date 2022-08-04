#!/bin/bash

# PHASE III Pt.1 - within gene directory

# create CSV files and plots; have plot making as an input option

# recieve genename as input: csv_plots.sh -a <gene_name> -b
while getopts "a:c" option
do 
	case $option in
		a) gene=$OPTARG;;
		c) plot_overlay=true;; #if user wants make overlay plot; optional
	esac
done

# get probe numbers
probes_list="$(ls $gene/*_final_pheno.txt)"

# For each probe, for each peak, run linear_to_csv.py <file.linear> <file.ld>
for PROBES in $probes_list
do
	# Prefix for downstream files: <gene>_probe_<ILMN_probe_number>
	prefix=${PROBES%.*}
	prefix=${prefix::${#prefix}-12}

	# univariate; make univariate files with: SNP, Position, P-value, Pk1_LD
	awk '{if ($2=="SNP") {print $2"\t"$3"\t"$9}}' "$prefix"_univariate.P5.assoc.linear > "$prefix"_univariate_tmp_linear_to_csv.txt 
	awk '{if ($5=="ADD") {print $2"\t"$3"\t"$9}}' "$prefix"_univariate.P5.assoc.linear >> "$prefix"_univariate_tmp_linear_to_csv.txt
	awk '{print $7}' "$prefix"_LD_peak1.ld > "$prefix"_temp_peak1.ld
	paste -d'\t' "$prefix"_univariate_tmp_linear_to_csv.txt "$prefix"_temp_peak1.ld > "$prefix"_univariate_linear_to_csv.txt
	cat "$prefix"_univariate_linear_to_csv.txt | tr -s '[:blank:]' ',' > "$prefix"_univariate.csv

	# peaks; loop through
	count=1 #peak number
	# For each peak w/in each probe, create csv files
	n="$(wc -l "$prefix"_ALPHA_SNPs_LIST.txt | awk '{print $1}')"

	while [ $count -le $n ]
	do
		# old: linear_to_csv.py "$prefix"_conditioned_peak"$count".P5.assoc.linear "$prefix"_LD_peak"$count".ld ; this script is no longer used.
		awk '{if ($2=="SNP") {print $2"\t"$3"\t"$9}}' "$prefix"_conditioned_peak"$count".P5.assoc.linear > "$prefix"_conditioned_peak"$count"_tmp_linear_to_csv.txt 
		awk '{if ($5=="ADD") {print $2"\t"$3"\t"$9}}' "$prefix"_conditioned_peak"$count".P5.assoc.linear >> "$prefix"_conditioned_peak"$count"_tmp_linear_to_csv.txt
		awk '{print $7}' "$prefix"_LD_peak"$count".ld > "$prefix"_temp_peak"$count".ld
		paste -d'\t' "$prefix"_conditioned_peak"$count"_tmp_linear_to_csv.txt "$prefix"_temp_peak"$count".ld > "$prefix"_conditioned_peak"$count"_linear_to_csv.txt
		cat "$prefix"_conditioned_peak"$count"_linear_to_csv.txt | tr -s '[:blank:]' ',' > "$prefix"_conditioned_peak"$count".csv
		count=$[$count+1]
	done

	# conditioned on all
	awk '{if ($2=="SNP") {print $2"\t"$3"\t"$9}}' "$prefix"_conditioned_all.P5.assoc.linear > "$prefix"_conditioned_all_tmp_linear_to_csv.txt 
	awk '{if ($5=="ADD") {print $2"\t"$3"\t"$9}}' "$prefix"_conditioned_all.P5.assoc.linear >> "$prefix"_conditioned_all_tmp_linear_to_csv.txt
	awk '{print $7}' "$prefix"_LD_peak1.ld > "$prefix"_temp_peak1.ld
	paste -d'\t' "$prefix"_conditioned_all_tmp_linear_to_csv.txt "$prefix"_temp_peak1.ld > "$prefix"_conditioned_all_linear_to_csv.txt
	cat "$prefix"_conditioned_all_linear_to_csv.txt | tr -s '[:blank:]' ',' > "$prefix"_conditioned_all.csv
done


rm $gene/*_linear_to_csv.txt

echo "CSV CONVERSION COMPLETE"

#---------- MAKE PLOTS -----------

# If user inputs plot option, then plot. If not, then do nothing.
if $plot_overlay; then
	for PROBES in $probes_list
	do
		# Prefix for downstream files: <gene>_probe_<ILMN_probe_number>
		prefix=${PROBES%.*}
		prefix=${prefix::${#prefix}-12}
		# (number files per probe) - 2 = # peaks
		n=$(ls $prefix*.csv | wc -l)
		peaks=$[$n-2]

		file_list=$(ls $prefix*.csv)

		if [ $peaks -le 10 ]
		then
			# input all files; parsed through via LZ_overlay_plots.R script; automatically makes plots for all
			Rscript LZ_overlay_plots.R $prefix $peaks "$prefix"_overlay "$prefix"_ALPHA_SNPs_LIST.txt $prefix*.csv

		else
			# if more than 10 peaks, tell user
			echo "Too many peaks to plot; more than 10 peaks"
		fi

	done
fi
