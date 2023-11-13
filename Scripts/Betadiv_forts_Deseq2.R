•	baseMean—The average of the normalized count values, dividing by size factors, taken over all samples.
•	log2FoldChange–The effect size estimate. This value indicates how much the gene or transcripts expression seems to have changed between the comparison and control groups. This value is reported on a logarithmic scale to base 2.
•	lfcSE–The standard error estimate for the log2 fold change estimate.
•	stat–The value of the test statistic for the gene or transcript.
•	pvalue–P-value of the test for the gene or transcript.
•	padj–Adjusted P-value for multiple testing for the gene or transcript.

# pakker jeg trenger
library(DESeq2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(apeglm)
library(DEGreport)
library(readxl)

# farger
color_palette_class_2 = c("#D1BBD7", "#AE76A3", "#882E72", "#1965B0",
                        "#5289C7", "#7BAFDE", "#4EB265", "#90C987",
                        "#CAE0AB", "#F7F056", "#F6C141", "#F1932D",
                        "#E8601C", "#DC050C", "#72190E")


# Laster inn filene mine og setter wd:

knitr::opts_knit$set(root.dir = "/Users/siljerosbak/Dropbox/Master - Silje/Data/") # hvor jeg henter dataene fra

metadata <- read_excel("Metadata1_master.xlsx") # laster inn metadataene
samples_df <- readRDS("samples_df.rds", refhook=NULL) # metadata

p_clean_gg_uten_2 <- readRDS("p_clean_gg_uten_2.rds", refhook=NULL) # phyloseq uten de negative kontroll filterene mine
p_clean_df <- psmelt(p_clean_gg_uten_2) # Lager data.frame



# Deseq2 på giftighet - Lager log-fold-change plot for å se forskjell i gruppene 
# Laster inn nødvendige pakker
dds <- phyloseq_to_deseq2(p_clean_gg_uten_2, ~ Giftighet) #bruker unrarefied fordi deseq2 gjør dette i neste steg + analysen blir gruppert etter giftighet
#dds$Giftighet

# keep <- rowSums(counts(dds) >=6) >=3 # must have 6 reads or more, in 3 samples or more
# table(keep)
# dds <- dds[keep,]

p_clean_gg_uten_2@sam_data %>%
  dplyr::group_by(Giftighet) %>%
  dplyr::summarise(n=n())

# Run the DESeq normalization and variance stabilizing transformation

gm_mean <- function(x, na.rm=TRUE){ 
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
geoMeans <- apply(counts(dds), 1, gm_mean)
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)

dds <- DESeq(dds, fitType="local") 

# assumptions
plotDispEsts(dds) #  check dispersion plot
resultsNames(dds) # checks for the groups produced. Using results function, contrast variable to specify which groups to compare.

#legge til grppe1, gruppe2 gruppe3 i metadata

#lager ulike lister som skal brukes til slutt:
results <- NULL
alge_data <- NULL
alge_frame <- NULL 

# gruppe 1 blir referanse, der gruppe to synes(hvor annerledes den er)
# alpha_gfitgihet <- 0.99
results_DEseq2_giftighet <- results(dds, contrast = c("Giftighet", "Potensielt_giftig", "Ikke-giftig")) # Get the differential abundance results


# Lage plot med log2 fold changes

alpha = 1 # adjust your p value here.
sigtab = results_DEseq2_giftighet[which(results_DEseq2_giftighet$padj <= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table((p_clean_gg_uten_2))[rownames(sigtab), ], "matrix")) # change ps_filt to your ps
head(sigtab) #checks
dim(sigtab) #checks
sigtab <- sigtab %>%
  rownames_to_column(var = "ASV")

# se på fordelinga for en spesifikt ASV
#plotCounts(dds, gene = "TTTGCCTACGGGTGGCTGCAGTGGGGAATCTTAGACAATGGGCGCAAGCCTGATCTAGCCATGCCGCGTGAGTGATGAAGGCCTTAGGGTCGTAAAGCTCTTTCGCCTGTGAAGATAATGACTGTAGCAGGTAAAGAAACCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGGGTTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGACATTTAAGTCAGAGGTGAAATCCCAGGGCTCAACCCTGGAACTGCCTTTGATACTGGGTGTCTTGAGTTCGAGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTGCGAAAGTGTGGGGAGCAAACAGGATTAGATACCCTAGTAGTC", intgroup = "Giftighet")

# asv_means<-data.frame()
# for(i in sigtab$ASV) {
#   Potensielt_giftig <- mean(counts(dds, normalized=TRUE)[i, which(dds$Giftighet == "Potensielt_giftig")])
#   Ikke-giftig <- mean(counts(dds, normalized=TRUE)[i, which(dds$Giftighet == "Ikke-giftig")])
#   asv = i
#   output = c("ASV", "Potensielt_giftig", "Ikke-giftig")
#   asv_means = rbind(asv_means, output)
# }
#--------- Plot
# Division and Order plot. Need change this manually if want to
# x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
# x = sort(x, TRUE)
# sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Order)) + # change here your level of plotting
  geom_point(size=3) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  scale_colour_manual(values=color_palette_class_2) + 
  geom_hline(yintercept = 0, linetype="solid") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x =element_text(size=10),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "right",
        strip.text.y = element_blank(),
        panel.spacing.y = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.border = element_blank(),
        strip.background = element_blank()
        # legend.position = "none"
        ) +
  ylab("log2 fold change - Ikke-giftig <---> Potensielt_giftig") +
  facet_grid(Order~., scales = "free", space= "free", switch = "y")
ggsave("deseq1.pdf", dpi=600)





# Gjør Deseq2 på fraksjon 
dds <- phyloseq_to_deseq2(p_clean_gg_uten_2, ~ Filter_strl_um) #bruker unrarefied fordi deseq2 gjør dette i neste steg + analysen blir gruppert etter giftighet

# Run the DESeq normalization and variance stabilizing transformation
gm_mean <- function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans <- apply(counts(dds), 1, gm_mean)
diagdds <- estimateSizeFactors(dds, geoMeans = geoMeans)

diagdds <- DESeq(diagdds, fitType="local") 


#fjerner listene som skal brukes på nytt
results <- NULL
alge_data <- NULL
alge_frame <- NULL 
sigtab <- NULL

results_DEseq2_fraksjon <- results(diagdds, contrast = c("Filter_strl_um", "Frittlevende(0,2um)", "Fastsittende(1um)")) # Get the differential


#Lage plot 
plotDispEsts(diagdds) #  check dispersion plot
resultsNames(diagdds) # checks for the groups produced. Using results function, contrast variable to specify which groups to compare.

alpha = 1 # adjust your p value here.
sigtab = results_DEseq2_fraksjon[which(results_DEseq2_fraksjon$padj <= alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table((p_clean_gg_uten_2))[rownames(sigtab), ], "matrix")) # change ps_filt to your ps

head(sigtab) #checks
dim(sigtab) #checks

# se på fordelinga for en spesifikt ASV
#plotCounts(dds, gene = "TGAACCTACGGGGGGCTGCAGTGAGGAATATTGGGCAATGGAGGCAACTCTGACCCAGCCATGCCGCGTGCAGGATGACGGCCCTATGGGTTGTAAACTGCTTTTATACAGGAATAAACCCCCGAACGTGTTCGGGGCTGAAGGTACTGTACGAATAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGGGTTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCGCGTAGGCGGACATTTAAGTCAGAGGTGAAATCCCAGGGCTCAACCCTGGAACTGCCTTTGATACTGGGTGTCTTGAGTTCGAGAGAGGTGAGTGGAATTCCGAGTGTAGAGGTGAAATTCGTAGATATTCGGAGGAACACCAGTGGCGAAGGCGGCTCACTGGCTCGATACTGACGCTGAGGTGCGAAAGTGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTC", intgroup = "Filter_strl_um")

sigtab <- sigtab %>%
  filter(log2FoldChange >= -10) # remove asv that only shows up once

#--------- Plot
# Division and Order plot. Need change this manually if want to
#x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
#x = sort(x, TRUE)
#sigtab$order = factor(as.character(sigtab$Order), levels=names(x))

ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Order)) + # change here your level of plotting
  geom_point(size=3) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+
  scale_colour_manual(values=color_palette_class_2) + 
  geom_hline(yintercept = 0, linetype="solid") +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.title.x =element_text(size=10),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "right",
        strip.text.y = element_blank(),
        panel.spacing.y = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.border = element_blank(),
        strip.background = element_blank()
        # legend.position = "none"
        ) +
  ylab("log2 fold change - Fastsittende(1um) <---> Frittlevende(0,2um") +
  facet_grid(Order~., scales = "free", space= "free", switch = "y")

ggsave("deseq2.pdf", dpi=600)
