#!/bin/bash

# PHASE II - within gene directory

# receeve genename as input: run_plink.sh -a <gene_name>
while getopts "a:" option
do 
	case $option in
		a) gene=$OPTARG;;
	esac
done

get_prefix="$(ls $gene/"$gene"*.fam)"
prefix=${get_prefix%.*}

plink=./plink

# Get all phenotype files; one per probe to iterate trhough
final_phenos="$(ls $gene/*_final_pheno.txt)"
echo $final_phenos

# do this protocol for each gene to run automatically through each probe
for PHENO_FILE in $final_phenos
do
	count=1 # peak number
	# Prefix for downstream files: <gene>_probe_<ILMN_probe_number>
	plink_prefix=${PHENO_FILE%.*}
	plink_prefix=${plink_prefix::${#plink_prefix}-12}
	# Find first peak
	$plink --bfile $prefix --linear --pheno $PHENO_FILE --all-pheno --allow-no-sex --out "$plink_prefix"_peak"$count"
	# Top peak is the smallest p-value; add top peak to ALPHA_SNPs_LIST
	top_peak="$(sort -k9 -g "$plink_prefix"_peak"$count".P5.assoc.linear | awk 'FNR == 2 {print $2}' > "$plink_prefix"_ALPHA_SNPs_LIST.txt)"
	top_pvalue="$(sort -k9 -g "$plink_prefix"_peak"$count".P5.assoc.linear | awk 'FNR == 2 {print $9}')"
	top_pvalue="$(echo "$top_pvalue" | tr '[a-z]' '[A-Z]')"
	top_pvalue="$(printf "%.14f\n" "$top_pvalue")"

	# Find all significant peaks, one at a time. Add each peak to ALPHA SNPs List
	if [[ $count -lt 20 ]]; then
	
		while [[ $(echo "$top_pvalue"'<'0.00005| bc -l) -eq 1 ]]
		do
			count=$[$count+1]
			$plink --bfile $prefix --linear --pheno $PHENO_FILE --all-pheno --condition-list "$plink_prefix"_ALPHA_SNPs_LIST.txt --allow-no-sex --out "$plink_prefix"_peak"$count"
			top_pvalue="$(awk '{ if ($5 == "ADD") {print}}' "$plink_prefix"_peak"$count".P5.assoc.linear | sort -k9 -g | awk '{ if ($9 != "NA") {print}}' | awk 'FNR == 2 {print $9}')"
			top_pvalue="$(echo "$top_pvalue" | tr '[a-z]' '[A-Z]')"
			
			# edit 5-09-19
        		top_pvalue="$(printf "%.14f\n" "$top_pvalue")"
	
			# if new peak is signigificant, add to SNP list
			if [[ $(echo "$top_pvalue"'<'0.00005| bc -l) -eq 1 ]]; then
				top_peak="$(awk '{ if ($5 == "ADD") {print}}' "$plink_prefix"_peak"$count".P5.assoc.linear | sort -k9 -g | awk '{ if ($9 != "NA") {print}}' | awk 'FNR == 1 {print $2}' >> "$plink_prefix"_ALPHA_SNPs_LIST.txt)"
			else
				rm "$plink_prefix"_peak"$count".P5.assoc.linear
			fi
			# Develop a checkpoint to select SNPs with p-values < 5.0 x 10^(-8) ; if peaks SNP is greater than or
			# equal to this p-value, stop. In other words, continue finding peaks and running plink until all sig. peaks
			# are found
			echo COUNT $count
			echo TOP $top_pvalue
		done
	fi
	# Run conditional analysis using all identified peaks. Need to have data/graphs for univariate, conditioned peaks1-n, and conditioned on all
	# 1. univariate; no conditioning
	cp "$plink_prefix"_peak1.P5.assoc.linear "$plink_prefix"_univariate.P5.assoc.linear
	# 2. condition on each peak individually
	# For each SNP in ALPHA_SNP_List:
	i=1
	for alpha_SNP in $(cat "$plink_prefix"_ALPHA_SNPs_LIST.txt)
	do
		echo $alpha_SNP
		# create a temp_condition_list; remove SNP of interest to condition on it
		sed '/'$alpha_SNP'/d' "$plink_prefix"_ALPHA_SNPs_LIST.txt > "$plink_prefix"_temp_"$i"_SNPs_LIST.txt
		# run plink2 with temp_condition_list
		$plink --bfile $prefix --linear --pheno $PHENO_FILE --all-pheno --condition-list "$plink_prefix"_temp_"$i"_SNPs_LIST.txt --allow-no-sex --out "$plink_prefix"_conditioned_peak"$i"
		
		# CHECKPOINT: document peaks that differ before/after conditional analysis
		beta_SNP="$(awk '{ if ($5 == "ADD") {print}}' "$plink_prefix"_conditioned_peak"$i".P5.assoc.linear | awk '{ if ($9 != "NA") {print}}' | sort -k9 -g | awk 'FNR == 1 {print $2}' | head)"
		if [ "$alpha_SNP" != "$beta_SNP" ]
		then
			echo "$plink_prefix" Peak"$i" Alpha SNP: "$alpha_SNP" Beta SNP: "$beta_SNP" >> "$plink_prefix"_DISAGREED_peak"$i".txt
			echo !!!SOME PEAKS DISAGREE!!! "$plink_prefix" Peak"$i" $alpha_SNP $beta_SNP
		fi

		# Calculate LD r^2 ; generate LD file for each top peak; PLINK.v.2.0 does not support this, only PLINK.v.1.9
		$plink --bfile $prefix --r2 --ld-snp "$alpha_SNP" --ld-window-kb 99999999 --ld-window 9999999 --ld-window-r2 0 --out "$plink_prefix"_LD_peak"$i"
		
		# CHECKPOINT: document peaks that differ before/after conditional analysis
		# next peak
		i=$[$i+1]
	done
	# 3. condition on all peaks
	$plink --bfile $prefix --linear --pheno $PHENO_FILE --all-pheno --condition-list "$plink_prefix"_ALPHA_SNPs_LIST.txt --allow-no-sex --out "$plink_prefix"_conditioned_all
done

rm $gene/*.logistic
rm $gene/*.log
rm $gene/*.nosex
#rm $gene/*.P1.assoc.linear
#rm $gene/*.P2.assoc.linear

echo DONE RUNNING PLINK
ls $gene/*DISAGREED*
