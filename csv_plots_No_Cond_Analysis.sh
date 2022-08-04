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

# put unconditioned csv files in a new dir to keep separate
mkdir not_conditioned/$gene
cp $gene/*_ALPHA_SNPs_LIST.txt not_conditioned/$gene

# get probe numbers
probes_list="$(ls $gene/*_final_pheno.txt)"

# For each probe, for each peak, combine .linear and .ld info
for PROBES in $probes_list
do
	# Prefix for downstream files: <gene>_probe_<ILMN_probe_number>
	prefix=${PROBES%.*}
	prefix=${prefix::${#prefix}-12}
	
	# univariate; make univariate files with: SNP, Position, P-value, Pk1_LD
	awk '{if ($2=="SNP") {print $2"\t"$3"\t"$9}}' "$prefix"_univariate.P5.assoc.linear > not_conditioned/"$prefix"_univariate_tmp_linear_to_csv.txt 
	awk '{if ($5=="ADD") {print $2"\t"$3"\t"$9}}' "$prefix"_univariate.P5.assoc.linear >> not_conditioned/"$prefix"_univariate_tmp_linear_to_csv.txt
	awk '{print $7}' "$prefix"_LD_peak1.ld > not_conditioned/"$prefix"_temp_peak1.ld
	paste -d'\t' not_conditioned/"$prefix"_univariate_tmp_linear_to_csv.txt not_conditioned/"$prefix"_temp_peak1.ld > not_conditioned/"$prefix"_univariate_linear_to_csv.txt
	cat not_conditioned/"$prefix"_univariate_linear_to_csv.txt | tr -s '[:blank:]' ',' > not_conditioned/"$prefix"_univariate.csv

	# peaks; loop through
	count=1 #peak number
	# For each peak w/in each probe, create csv files
	n="$(wc -l not_conditioned/"$prefix"_ALPHA_SNPs_LIST.txt | awk '{print $1}')"
	while [ $count -le $n ]
	do
		# old: linear_to_csv.py "$prefix"_peak"$count".P5.assoc.linear "$prefix"_LD_peak"$count".ld ; this script is no longer used.
		awk '{if ($2=="SNP") {print $2"\t"$3"\t"$9}}' "$prefix"_peak"$count".P5.assoc.linear > not_conditioned/"$prefix"_peak"$count"_tmp_linear_to_csv.txt 
		awk '{if ($5=="ADD") {print $2"\t"$3"\t"$9}}' "$prefix"_peak"$count".P5.assoc.linear >> not_conditioned/"$prefix"_peak"$count"_tmp_linear_to_csv.txt
		awk '{print $7}' "$prefix"_LD_peak"$count".ld > not_conditioned/"$prefix"_temp_peak"$count".ld
		paste -d'\t' not_conditioned/"$prefix"_peak"$count"_tmp_linear_to_csv.txt not_conditioned/"$prefix"_temp_peak"$count".ld > not_conditioned/"$prefix"_peak"$count"_linear_to_csv.txt
		cat not_conditioned/"$prefix"_peak"$count"_linear_to_csv.txt | tr -s '[:blank:]' ',' > not_conditioned/"$prefix"_peak"$count".csv
		count=$[$count+1]
	done

	# conditioned on all
	awk '{if ($2=="SNP") {print $2"\t"$3"\t"$9}}' "$prefix"_conditioned_all.P5.assoc.linear > not_conditioned/"$prefix"_conditioned_all_tmp_linear_to_csv.txt 
	awk '{if ($5=="ADD") {print $2"\t"$3"\t"$9}}' "$prefix"_conditioned_all.P5.assoc.linear >> not_conditioned/"$prefix"_conditioned_all_tmp_linear_to_csv.txt
	awk '{print $7}' "$prefix"_LD_peak1.ld > not_conditioned/"$prefix"_temp_peak1.ld
	paste -d'\t' not_conditioned/"$prefix"_conditioned_all_tmp_linear_to_csv.txt not_conditioned/"$prefix"_temp_peak1.ld > not_conditioned/"$prefix"_conditioned_all_linear_to_csv.txt
	cat not_conditioned/"$prefix"_conditioned_all_linear_to_csv.txt | tr -s '[:blank:]' ',' > not_conditioned/"$prefix"_conditioned_all.csv
done


rm not_conditioned/$gene/*_linear_to_csv.txt

echo "CSV CONVERSION COMPLETE"

#---------- MAKE PLOTS -----------

# If user inputs plot option, then plot. If not, then do nothing.
if [ $plot ]; then

	echo "Creating Locus Zoom Plots..."
	# plot every csv file for each probe. Put into new directory: gene_probe_lz_plots

	uni_list="$(ls not_conditioned/$gene/*univariate*.csv)"
	for UNI in $uni_list
	do
		out=$UNI
		out=${out::${#out}-4}
		echo $out
		LZ_plot_R.R $UNI $out
	done

	peak_list="$(ls not_conditioned/$gene/*peak*.csv)"
	for PEAK in $peak_list
	do
		out=$PEAK
		out=${out::${#out}-4}
		peak_num=${out: -1}
		probe_num=${out::${#out}-6}
		echo PROBE NUM $probe_num
		# get peak SNP
		peak_snp="$(awk 'NR=='$peak_num'{print$1}' "$probe_num"_ALPHA_SNPs_LIST.txt)"
		echo TESTING "$probe_num"_ALPHA_SNPs_LIST.txt
		LZ_plot_R.R $PEAK $out $peak_snp
	done

	cond_list="$(ls not_conditioned/$gene/*conditioned_all*.csv)"
	for CONDITIONED in $cond_list
	do
		out=$CONDITIONED
		out=${out::${#out}-4}
		echo $out
		LZ_plot_R.R $CONDITIONED $out
	done

echo "PLOTS MADE"	
fi

# If user wants to make overlaying plots:
if $plot_overlay; then
	for PROBES in $probes_list
	do
		# Prefix for downstream files: <gene>_probe_<ILMN_probe_number>
		prefix=${PROBES%.*}
		prefix=${prefix::${#prefix}-12}
		echo PREFIX $prefix
		# (number files per probe) - 2 = # peaks
		n=$(ls not_conditioned/$prefix*.csv | wc -l)
		peaks=$[$n-2]

		file_list=$(ls not_conditioned/$prefix*.csv)

		if [ $peaks -le 9 ]
		then
			# input all files; parsed through via LZ_overlay_plots.R script; automatically makes plots for all
			Rscript ./LZ_overlay_plots.R not_conditioned/$prefix $peaks not_conditioned/"$prefix"_overlay not_conditioned/"$prefix"_ALPHA_SNPs_LIST.txt not_conditioned/$prefix*.csv

		else
			# if more than 9 peaks, tell user
			echo "Too many peaks to plot; more than 9 peaks"
		fi

	done
fi
