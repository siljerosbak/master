# Lage phyloseq objekt og rarefisere data


# Pakker jeg trenger
library(phyloseq)
library(ggplot2)      # graphics
library(readxl)       # necessary to import the data from Excel file
library(dplyr)        # filter and reformat data frames
library(tibble)       # Needed for converting column to row names
library(knitr)
library(tidyverse)
library(micorbiome)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")



# Laster inn filene mine og setter wd:

knitr::opts_knit$set(root.dir = "/Users/siljerosbak/Dropbox/Master - Silje/Data/") # hvor jeg henter dataene fra
file.path("Silje - Master/Data", "seqtab.rds")


seqtab <- readRDS("seqtab.rds") # laster inn sekvensfil
tax_final <- readRDS("tax_final.rds") # laster inn taksonomi-fil
metadata <- read_excel("Metadata1_master.xlsx") # laster inn metadataene


# Lager phylo-object:

y <- readRDS("seqtab.rds", refhook=NULL)
otu <- t(y) # transformerer slik at den får sekvensene som ege kolonne/rad?

x <- t(seqtab)
tax <- readRDS("tax_final.rds", refhook=NULL)

write_xlsx(tox, "/Users/siljerosbak/Dropbox/Master - Silje/Data/tox_asv.xlsx")

otu <- otu_table(x, taxa_are_rows = T) # kaller seqtab-filen (x) for otu
tax <- tax_table(as.matrix(tax_final)) # kaller tax-final-filen for tax 

otu_mat <- as.matrix(otu) 

samples_df <- metadata %>% # kaller metadataene for samples_df
  tibble::column_to_rownames("Prøve")

samples_df$Sample <- row.names(samples_df) # legger til en kolonne der Samples blir en egen 

OTU <- otu_table(otu, taxa_are_rows = T) 
otu_filt <- OTU[rowSums(OTU) !=0,] #velger ut de radene der summen ikke er lik 0, dermed beholder vi de som ikke er 0
TAX <- tax_table(tax)

samples <- sample_data(samples_df)

phylo <- phyloseq(otu, tax, samples) # slår de to tabellene sammen i phyloseq
phylo_clean <- subset_taxa(phylo, Kingdom  %in% c("Bacteria","Archea")) # jeg vil se på bakterier og arker, men er ingen arker 
phylo_clean <- subset_taxa(phylo_clean, Order != "Chloroplast") # clean for chloroplasts 
phylo_clean <- subset_taxa(phylo_clean, Family != "Mitochondria") # er ingen mitokondrier 
phylomelt <- psmelt(phylo_clean) # lager et stort datasett ved å slå sammen otu, tax og metadata 





# hvor mange 
# Antall reads i phylo
rowSums(otu_table(subset_taxa(phylo))) %>% as.tibble()

phylo_abund <- psmelt(otu_table(subset_taxa(phylo))) %>% as.tibble()

sum(phylo_abund$Abundance)


# Antall reads som er Chloroplast
rowSums(otu_table(subset_taxa(phylo, Order = "Chloroplast"))) %>% as.tibble()

chlor_abund <- psmelt(otu_table(subset_taxa(phylo, Order == "Chloroplast"))) %>% as.tibble()

sum(chlor_abund$Abundance)


# Antall reads i phylo_clean
rowSums(otu_table(subset_taxa(phylo_clean))) %>% as.tibble()

phylo_clean_abund <- psmelt(otu_table(subset_taxa(phylo_clean))) %>% as.tibble()

sum(phylo_clean_abund$Abundance)




# Lager rarefraction curve

otu_clean <- otu_table(phylo_clean, taxa_are_rows = TRUE)
class(otu_clean) <- "matrix"
vegan::rarecurve(t(otu_clean),
                 sample = 5000,
                 step = 20)

# legger til disse linjene for å tilpasse figuren
Rare_tidy <- vegan::rarecurve(t(otu_clean),
                              tidy = TRUE)


ggplot(data = Rare_tidy, aes(x = Sample, y = Species, group = Site, col = Site)) + 
  geom_line(lty = 2) + 
  geom_vline(xintercept = 5000, lty = 1, lwd = 0.3) + 
  theme(legend.position = "none") +
  labs(x = "Antall avlesninger", y = "Antall ASVer") +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=8),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        panel.border = element_rect(color = "gray", fill = NA, size = .5)) 
#ggsave("rare_plot.pdf", dpi=600)



# fjerner de to negative kontrollene: 

p_clean_gg_uten_2 <- subset_samples(phylo_clean, !(Art %in% c("NA"))) # Fjerner negative 
df_phylo <- psmelt(p_clean_gg_uten_2)



# Rarefying samples:
# plukker ut de prøvene som har mindre enn 5000 reads
stor_ps <- subset_samples(p_clean_gg_uten_2, Sample !="SR-24" & Sample !="SR-27" & Sample !="SR-29" & Sample !="SR-37")

liten_ps <- subset_samples(p_clean_gg_uten_2, Sample %in% c("SR-24","SR-27","SR-29","SR-37"))

# rarifisere
phylo_5000 <- prune_samples(sample_sums(stor_ps)>= 5000, stor_ps)

set.seed(1) # må det være et tall?
rarefied_ps_1 <- rarefy_even_depth(stor_ps, sample.size = min(sample_sums(phylo_5000)),
                                   rngseed = TRUE, replace = TRUE, trimOTUs = FALSE, verbose = TRUE)  #Normalizing species

rarefied_ps <- merge_phyloseq(rarefied_ps_1, liten_ps)

p_rarefied_df <- psmelt(rarefied_ps) # lager en dataframe

length(unique(p_rarefied_df$OTU)) # lengden på antall ASVer 

#summarize_phyloseq(rarefied_ps) # for å sjekke at det er like mange reads i alle prøvene mine



# Hvor mange ASVer/reads i hhv rarified of ikke-rarified?

#print(phylomelt) # datasettet før sjeldenhetsmpling og subsamling 
length(unique(phylomelt$OTU))  #7241


# Antall reads i phylo_clean_df
sum(phylomelt$Abundance)



#print(p_rarefied) - rarified datasett
length(unique(p_rarefied_df$OTU)) #7241


# Antall reads i phylo_clean_df
sum(p_rarefied_df$Abundance)




# Lager RDS objekt til senere bruk

saveRDS(samples_df, file = "samples_df.rds") #metadata
saveRDS(rarefied_ps, file = "rarefied_ps.rds") #phyloseq som er rarifisert
saveRDS(p_rarefied_df, file = "p_rarefied_df.rds") # phyloseq som er blit data frame
saveRDS(p_clean_gg_uten_2, file = "p_clean_gg_uten_2.rds") # phyloseq uten de negative kontroll filterene mine


