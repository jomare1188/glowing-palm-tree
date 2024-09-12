### REFORMAT MODULES OUTPUT
count=1
file=../results/networks/mcl_cane_p0.6_i2.7.mcl
rm ${file}.formated.csv

for line in $(cat $file)
        do
                sed -n "${count}p" $file | awk -F" " -v c="$count" '{for(i=1; i<=NF; i++) {print $i,c}}' >> ${file}.formated.csv
                let count++
        done

