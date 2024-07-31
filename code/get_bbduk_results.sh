# Script to collect bbduk results from multiple accessions and make a unified table
# format table accession, genotype, total, matched, %matched.
rm -rf pretty_sorghum_table_bbduk.csv

for genotype in $(cat accessions/genotypes.txt)
do
	results_dir=../data/${genotype}/2_trimmed
	for sample in $(sort -u ids_sra/${genotype}_ids.txt)
	do
		total=$(grep -w "Total" ${results_dir}/${sample}_trimmed.stats | cut -f2)
		matched=$(grep -w "Matched" ${results_dir}/${sample}_trimmed.stats | cut -f2)
		percentage=$(grep -w "Matched" ${results_dir}/${sample}_trimmed.stats | cut -f3 | sed s/%//)
		echo ${sample},${genotype},${total},${matched},${percentage} >> pretty_sorghum_table_bbduk.csv
	done
done
