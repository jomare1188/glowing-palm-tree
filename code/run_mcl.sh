in_sorghum=../results/networks/triplets_p60_sorghum_fancy.tsv
out_sorghum=../results/networks/mcl_p0.6_i2.7.mcl

in_cane=../results/networks/triplets_p60_cane_fancy.tsv
out_cane=../results/networks/mcl_cane_p0.6_i2.7.mcl

#mcl $in_sorghum -I 2.7 -te 50 --abc -o ${out_sorghum}
mcl $in_cane -I 2.7 -te 50 --abc -o ${out_cane}
