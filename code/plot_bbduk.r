raw <- read.table("pretty_table_bbduk.csv", col.names = c("Accession", "Genotype", "Total", "Matched", "Percentage"), header = F, sep = ",")
library(dplyr)
library(ggplot2)
library(svglite)
library(ggrepel)

raw <- raw[complete.cases(raw),]
find_outlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

df <- raw %>%
        group_by(Genotype) %>%
        mutate(outlier = ifelse(find_outlier(Percentage), Accession, NA))


ggplot(data = df, aes (x = Genotype, y = Percentage))+
  geom_boxplot()+
# geom_text_repel(aes(label=outlier), na.rm=TRUE) +
  ggtitle("Sugarcane")+
  xlab("Genotype") + 
  ylab("% Removed Sequences") +
  ylim(0, 100) +
  theme_bw() +
  theme( panel.grid.major = element_blank(),
         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 14, face="bold"),
         axis.text.y = element_text(size =14, face="bold"),
         plot.title = element_text(size=22), 
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.line = element_line(colour = "black")
  )

ggsave(filename = "bbduk_removed_seqs_sugarcane.svg", device = "svg", units = "in", width = 6.05*1.2, height = 3.75*1.2)
ggsave(filename = "bbduk_removed_seqs_sugarcane.png", device = "png", units = "in", width = 6.05*1.2, height = 3.75*1.2)


#make big pretty table with all preprocessing
sugarcane <- read.table("../results/bbduk/pretty_table_bbduk.csv", col.names = c("Accession", "Genotype", "Total", "bbduk_removed", "bbduk_percentage_removed"), header = F, sep = ",")
sorghum <- read.table("../results/bbduk/pretty_sorghum_table_bbduk.csv", col.names = c("Accession", "Genotype", "Total", "bbduk_removed", "bbduk_percentage_removed"), header = F, sep = ",")

sugarcane$Species <- "Sugarcane"
sorghum$Species <- "Sorghum"

sugarcane <- sugarcane %>% mutate(bbduk_remaining = Total - bbduk_removed)
sorghum <- sorghum %>% mutate(bbduk_remaining = Total - bbduk_removed)

salmon_sugarcane <- read.table("../results/quant/pretty_table_sugarcane.csv", col.names = c("Accession", "Genotype", "salmon_mapped_percentage_respect_bbduk"), header = F, sep = ",")
salmon_sorghum <- read.table("../results/quant/pretty_table_sorghum.csv", col.names = c("Accession", "Genotype", "salmon_mapped_percentage_respect_bbduk"), header = F, sep = ",")

salmon_sugarcane$Species <- "Sugarcane"
salmon_sorghum$Species <- "Sorghum"

merged_sugarcane <- merge(sugarcane, salmon_sugarcane, by = c("Accession", "Genotype", "Species"))
merged_sorghum <- merge(sorghum, salmon_sorghum, by = c("Accession", "Genotype", "Species"))


merged_sugarcane <- merged_sugarcane %>% mutate(salmon_remaining = round(bbduk_remaining*salmon_mapped_percentage_respect_bbduk/100)) %>% mutate(final_percentage_remaining_respect_total = round(100*(salmon_remaining/Total),4))
merged_sorghum <- merged_sorghum %>% mutate(salmon_remaining = (bbduk_remaining*salmon_mapped_percentage_respect_bbduk)/100) %>% mutate(final_percentage_remaining_respect_total = 100*(salmon_remaining/Total))

final <- rbind(merged_sugarcane, merged_sorghum)

write.table(final, "../results/summary_cleaning.csv", quote = F, row.names = F, col.names = T, sep = ",")


