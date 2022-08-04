#!/bin/bash

# receive genename as input: extract_phenotypes.sh -a <gene_name>
while getopts "a:" option
do 
	case $option in
		a) gene=$OPTARG;;
	esac
done

#Create OG Phenotype Files for Each Probe per Gene

# count how many times gene appears in phenotype file; make that many files
n="$(awk -v a="$gene" '{if ($2 == a) {print}}' phenotype_for_chosen_genes | wc -l)"
#echo $n

# make directory for gene 
mkdir $gene

# Move Genotype Files to New $gene Directory 
cp ./genotype/$gene"_"* ./$gene

# generate correct phenotype format files for PLINK
i=1
while [ $i -le $n ]
do
	probe="$(awk -v gene="$gene" '{if ($2==gene) {print}}' phenotype_for_chosen_genes | awk -v n="$i" '{if (NR==n) {print $3}}')"
	echo $probe
	awk -v p=$probe '{if ($3 == p) {print}}' phenotype_for_chosen_genes | cut -f 4- > temp_exp_data.txt
	./transpose.sh -a temp_exp_data.txt -b trans_temp_exp_data.txt
	paste -d"\t" ./$gene/$gene"_"*fam trans_temp_exp_data.txt > ./$gene/"$gene""_probe_""$probe"_final_pheno.txt
	i=$[$i+1]
done
