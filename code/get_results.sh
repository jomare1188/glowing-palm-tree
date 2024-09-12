# all_cds, all_transcripts, longest_cds_per_orthogroup, ref
# format table accession,reference,mapping_rate
rm -rf pretty_table.csv

for genotype in $(cat accessions/genotypes.txt)
do
	results_dir=../data/${genotype}/3_salmon_quant_longest
	for sample in $(sort -u ids_sra/${genotype}_ids.txt)
	do
		mapping_rate=$(grep -w "Mapping rate" ${results_dir}/${sample}/logs/salmon_quant.log | rev | cut -f1 -d" " | rev | sed "s/%//g")
		echo ${sample},${genotype},${mapping_rate} >> pretty_table_sorghum.csv
	done
done
