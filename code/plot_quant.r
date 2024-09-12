raw <- read.table("pretty_table_cane.csv", col.names = c("Accession", "Genotype", "Mapping"), header = F, sep = ",")
library(dplyr)
library(ggplot2)

ggplot(data = raw, aes (x = Genotype, y = Mapping))+
geom_boxplot()+
  ggtitle("Sugarcane")+
  xlab("Genotype") + 
  ylab("Mapping %") +
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

ggsave(filename = "quant_cane.svg", device = "svg", units = "in", width = 6.05*1.2, height = 3.75*1.2)
ggsave(filename = "quant_cane.png", device = "png", units = "in", width = 6.05*1.2, height = 3.75*1.2)

