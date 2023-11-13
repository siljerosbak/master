# Taksonomi 

# pakker jeg trenger
library(ggplot2)
library(ggVennDiagram)
library(ggvenn)
library(grid)
library(writexl)
library(treemap)
library(readxl)
library(scales)

# Laster inn filene mine og setter wd:

knitr::opts_knit$set(root.dir = "/Users/siljerosbak/Dropbox/Master - Silje/Data/") # hvor jeg henter dataene fra

# Laster inn phyloseq objektene jeg trenger

metadata <- read_excel("Metadata1_master.xlsx") # laster inn metadataene
samples_df <- readRDS("samples_df.rds", refhook=NULL) # metadata
rarefied_ps <- readRDS("rarefied_ps.rds", refhook=NULL) # rarefied phyloseq (2028)
p_rarefied_df <- readRDS("p_rarefied_df.rds", refhook=NULL) # rarefied data frame



# Taksonomisk identifikasjon
## sjekke hvor mange fyla det er i datasettet 
antall_fyla <- p_rarefied_df$Phylum
length(antall_fyla)
length(unique(antall_fyla)) # 3
list(unique(antall_fyla)) # "Proteobacteria"   "Bacteroidota"     "Actinobacteriota"

## sjekke hvor mange ordener det er i datasettet 
antall_ordener <- p_rarefied_df$Order
length(antall_ordener)
length(unique(antall_ordener)) # 14 - som stemmer overens med figuren

## sjekke hvor mange fyla det er i datasettet 
antall_genus <- p_rarefied_df$Genus
length(antall_genus)
length(unique(antall_genus)) # 54, men en er NA - så det er 53
list(unique(antall_genus))

#print(p_rarefied_df$Phylum)



# legge inn farger og tema
color_palette_class_2 = c("#D1BBD7", "#AE76A3", "#882E72", "#1965B0",
                          "#5289C7", "#7BAFDE", "#4EB265", "#90C987",
                          "#CAE0AB", "#F7F056", "#F6C141", "#F1932D",
                          "#E8601C", "#DC050C", "#72190E")


mitt_tema <- theme(axis.text = element_text(size=10),
                   axis.title = element_text(size=14, face="bold"),
                   legend.text = element_text(size=12),
                   plot.background = element_rect(fill = "white"),
                   panel.background = element_rect(fill = "white"),
                   axis.line.x = element_line(color = "grey")) 



class_colors_2 <- setNames(color_palette_class_2, levels(p_rarefied_df$Order))


# Create the treemap with the transparent colors
treemap(p_rarefied_df,
        index=c("Phylum", "Order"),
        vSize="Abundance",
        type = "index",
        palette = color_palette_class_2)

sequences_df <- p_rarefied_df %>% # Lager d.f med de unike sekvensene
  group_by(Phylum, Order) %>%
  summarise(unique_sequences = n_distinct(OTU))

treemap(sequences_df, 
        index=c("Phylum","Order"),
        vSize="unique_sequences",
        vColor=y,
        type="index",
        palette = color_palette_class_2)
#vp = vplayout(1,2))
#Endrer farger i Inkskape


# lager bar plot med rarified for en oversikt over taksonomien - de tre rekkene
ggplot(p_rarefied_df, aes(fill=Phylum,y=Abundance, x = as.factor(Label))) + 
  geom_bar(stat = "identity", position = "fill", alpha = 0.6) +
  facet_grid(Filter_strl_um ~ Giftighet,
             scales = c("free"),
             space = "free_y",
             drop = TRUE) +
  coord_flip() +
  xlab("Prøver") +
  ylab("Andel avlesninger [%]") +
  mitt_tema +
  scale_fill_manual(values =  c("#DC050C", "#F7F056", "#F1932D")) +
  scale_y_continuous(labels = scales::percent_format(scale = 100))  # Format y-axis as percentages
#ggsave("ggplot2.png", dpi=600)




# lager bar plot med rarified for en oversikt over taksonomien - orden
ggplot(p_rarefied_df, aes(fill=Order, y=Abundance, x=as.factor(Label))) + 
  geom_bar(stat = "identity", position = "fill", alpha = 0.6) +
  facet_grid(Filter_strl_um ~ Giftighet,
             scales = "free",  # No need to use a vector here
             space = "free_y",
             drop = TRUE) +
  coord_flip() +
  xlab("Prøver") +
  ylab("Andel avlesninger [%]") +
  mitt_tema + 
  scale_fill_manual(values = class_colors_2) +
  scale_y_continuous(labels = scales::percent_format(scale = 100))  # Format y-axis as percentages
#ggsave("ggplot1.pdf", dpi=600)




# Sort the data by Frequency and select the top 10 values
Genus_unique <- p_rarefied_df %>%
  distinct(Genus, .keep_all = TRUE)

Genus_unique_2 <- phyloseq::tax_glom(rarefied_ps, "Genus")

Genus_unique_3 <-  (Genus_unique_2_df %>% select(Genus, Sample, Abundance))

Genus_unique_4 <- aggregate(Abundance ~ Genus, data = Genus_unique_3, FUN=sum) # tabell i appendiks

#write_xlsx(Genus_unique_4, "/Users/siljerosbak/Dropbox/Master - Silje/Data/Genus_unique_4.xlsx")



# lage et søylediagram som viser de 15 mest abundante slektene - %-vis

taxonomy = "Genus"

classGlom = phyloseq::tax_glom(rarefied_ps, taxrank = taxonomy)

taxon_table = phyloseq::otu_table(classGlom)
tax_matrix = as(phyloseq::tax_table(classGlom), 'matrix')
rownames(taxon_table) = tax_matrix[,taxonomy]
tax_table = prop.table(taxon_table, margin = 2)*100
tax_table <- tax_table[order(rowSums(-taxon_table)),]
rowSums(taxon_table)
dim(taxon_table)
tax_table <- tax_table[1:15,]
rownames(tax_table)

dom1 <- apply(taxon_table,1,which.max)
list1 <- rownames(as.data.frame(taxon_table))[dom1]
table(list1)
# calculate proportions based on phylum data

melted_tax_table = reshape2::melt(tax_table)
melted_tax_table <- arrange(melted_tax_table, Var1, desc(value))
melted_tax_table$Var2 <- factor(melted_tax_table$Var2 , levels = unique(melted_tax_table$Var2))
metadata = phyloseq::sample_data(classGlom)
treatment_id <- metadata$Label[match(melted_tax_table$Var2, rownames(metadata))]

class_colors <- setNames(color_palette_class_2, levels(melted_tax_table$Var1))

#cairo_ps(file.path(figs_dir, "genus_plot.eps"),
#width=120/25.4,
#height=160/25.4,
#pointsize = 4,
#bg = FALSE,
#fallback_resolution = 300)
ggplot(melted_tax_table, aes(x = treatment_id, y = value, fill = Var1))+
  geom_bar(stat = "identity", position = "stack", alpha = .5) +
  guides(fill = guide_legend(title = taxonomy)) +
  coord_flip() +
  mitt_tema +
  xlab("Prøver") +
  ylab("Andel avlesninger [%]") +
  scale_fill_manual(values = class_colors) +
  scale_x_discrete(limits = rev(levels(treatment_id)))





# Venndiagram - hvor mange bakterier finnes kun i de giftige, som ikke finnes i de ugiftige
df_venn <- p_rarefied_df %>% filter(Abundance > 0)

venn_giftighet <- split(df_venn$OTU, df_venn$Giftighet)
ggVennDiagram(venn_giftighet, label_alpha = 0) +
  scale_fill_gradient(low ="#CAE0AB", high = "#E8601C") +
  labs(fill = "Andel ASVer") # venndiagram over ASV
ggsave("venn1.pdf", dpi=600)


venn_genus <- split(df_venn$Genus, df_venn$Giftighet)
ggVennDiagram(venn_genus, label_alpha = 0)+
  scale_fill_gradient(low ="#CAE0AB", high = "#E8601C") +
  labs(fill = "Andel slekter") # venndigram over slekt
ggsave("venn2.pdf", dpi=600)

venn_art <- split(df_venn$OTU, df_venn$Art)
ggVennDiagram(venn_art, label_alpha = 0)+
  scale_fill_gradient(low ="#CAE0AB", high = "#E8601C") +
  labs(fill = "Andel ASV") # venndigram over de seks ulike artene
ggsave("venn3.pdf", dpi=600)

df_chryso <- p_rarefied_df %>% filter(Abundance > 0,
                                      Art == "Chrysochromulina leadbeateri")
venn_chryso <- split(df_chryso$OTU, df_chryso$Stamme...4)
ggVennDiagram(venn_chryso, label_alpha = 0) +
  scale_fill_gradient(low ="#CAE0AB", high = "#E8601C") +
  labs(fill = "Andel ASVer") # venndiagram over ASV
ggsave("venn4.pdf", dpi=600)

venn_chryso_genus <- split(df_chryso$Genus, df_chryso$Stamme...4)
ggVennDiagram(venn_chryso_genus, label_alpha = 0)+
  scale_fill_gradient(low ="#CAE0AB", high = "#E8601C") +
  labs(fill = "Andel slekter") # venndiram over slekt
ggsave("venn5.pdf", dpi=600)


## Fordelingen av andel avlesninger
venn_abundans <- split(df_venn$OTU, df_venn$Giftighet)
lapply(venn_abundans, FUN=sum) # for å se hvor mange reads det er i de to gruppene

((64331 * 100)/(64331 + 128963)) # 33%

((128963 * 100)/(64331 + 128963)) # 66%

#bare for å sjekke:
64331 + 128963  # Tot reads = 193294 



## For å finne de unike slekten som kun finnes i de *giftige* algene: 
tox <- filter(p_rarefied_df, Giftighet == "Potensielt_giftig", Abundance >= 1)
non_tox <- filter(p_rarefied_df, Giftighet == "Ikke-giftig", Abundance >= 1)

venn <- list(TOX = tox$OTU, NON_TOX = non_tox$OTU)
venn <- list(TOX = tox$OTU, NON_TOX = non_tox$OTU)
summary(p_rarefied_df$Abundance) # en liste som viser at flesteparten av prøvene mine har observasjonen 0 - derfor fjerner jeg dette fra datasettet ved å sette "Abundance >=1" i linja over

venn_2 <- list(TOX = tox$Genus, NON_TOX = non_tox$Genus)
venn_3 <- list(TOX = tox$Order, NON_TOX = non_tox$Order)


groupA_unique <- setdiff(venn_2$TOX, venn_2$NON_TOX)

groupA_df <- subset(p_rarefied_df, Genus %in% groupA_unique, select = c("Genus", "Abundance"))

unique_tox_df <- groupA_df %>%
  distinct(Genus, .keep_all = TRUE)

# lage tabell til excel:
write_xlsx(unique_tox_df, "/Users/siljerosbak/Dropbox/Master - Silje/Data/unike_slekter_giftige.xlsx")



# For å finne de unike slekten som kun finnes i de *ikke-giftige* algene: 
groupB_unique <- setdiff(venn_2$NON_TOX, venn_2$TOX)

groupB_df <- subset(p_rarefied_df, Genus %in% groupB_unique, select = c("Genus", "Abundance"))

unique_nontox_df <- groupB_df %>%
  distinct(Genus, .keep_all = TRUE)
distinct(Genus, .keep_all = TRUE)
# lage tabell til excel:
#write_xlsx(unique_nontox_df, "/Users/siljerosbak/Dropbox/Master - Silje/Data/unike_slekter_ikkegiftige.xlsx")


# For å finne de slektene som  finnes i *BEGGE* giftighetene: 

groupC_unique <- intersect(venn_2$NON_TOX, venn_2$TOX)

groupC_df <- subset(p_rarefied_df, Genus %in% groupC_unique, select = c("Genus", "Abundance"))

begge_grupper_df <- groupC_df %>%
  distinct(Genus, .keep_all = TRUE)
distinct(Genus, .keep_all = TRUE)
# lage tabell til excel:
#write_xlsx(begge_grupper_df, "/Users/siljerosbak/Dropbox/Master - Silje/Data/begge_giftigheter.xlsx")






# Venndiagram - hvor mange bakterier finnes på algecellene vs fritt i kultur
På_algecellen <- filter(p_rarefied_df, Filter_strl_um == "Fastsittende(1um)", Abundance >= 1)
Fritt_i_kultur <- filter(p_rarefied_df, Filter_strl_um == "Frittlevende(0,2um)", Abundance >= 1)


ggVennDiagram(list(Fastsittende = På_algecellen$OTU, Frittlevende = Fritt_i_kultur$OTU), 
              label_alpha = 0) + 
  scale_fill_gradient(low ="#CAE0AB", high = "#E8601C") + 
  labs(fill = "Andel slekter")
ggsave("ASV_filter_venn.pdf", dpi=600)


ggVennDiagram(list(Fastsittende = På_algecellen$Genus, Frittlevende = Fritt_i_kultur$Genus), 
              label_alpha = 0) + 
  scale_fill_gradient(low ="#CAE0AB", high = "#E8601C") +
  labs(fill = "Andel slekter")
ggsave("ASV_filter_venn2.pdf", dpi=600)




## For å finne de unike slekten som kun finnes *på algene*: 
alge_venn <- list(ALGE = På_algecellen$Genus, KULTUR = Fritt_i_kultur$Genus)

alge_unique <- setdiff(alge_venn$ALGE, alge_venn$KULTUR)

alge_uniqe_df <- subset(p_rarefied_df, Genus %in% alge_unique, select = c("Genus", "Abundance"))

unike_slekter_fastsittende <- alge_uniqe_df %>%
  distinct(Genus, .keep_all = TRUE)
# lage tabell via excel:
#write_xlsx(unike_slekter_fastsittende, "/Users/siljerosbak/Dropbox/Master - Silje/Data/unike_slekter_fastsittende.xlsx")




## For å finne de unike slekten som kun finnes *fritt i kultur*: 
alge_venn <- list(ALGE = På_algecellen$Genus, KULTUR = Fritt_i_kultur$Genus)

alge_unique_2 <- setdiff(alge_venn$KULTUR, alge_venn$ALGE)

alge_uniqe_df_2 <- subset(p_rarefied_df, Genus %in% alge_unique_2, select = c("Genus", "Abundance"))

unike_slekter_frittlevende <- alge_uniqe_df_2 %>%
  distinct(Genus, .keep_all = TRUE)
# lage tabell via excel:
#write_xlsx(unike_slekter_frittlevende, "/Users/siljerosbak/Dropbox/Master - Silje/Data/unike_slekter_frittlevende.xlsx")




## For å finne de slektene som  finnes i *BEGGE*: 
alge_venn <- list(ALGE = På_algecellen$Genus, KULTUR = Fritt_i_kultur$Genus)

alge_unique_3 <- intersect(alge_venn$KULTUR, alge_venn$ALGE)

alge_uniqe_df_3 <- subset(p_rarefied_df, Genus %in% alge_unique_3, select = c("Genus", "Abundance"))

felles_slekter_fraksjon <- alge_uniqe_df_3 %>%
  distinct(Genus, .keep_all = TRUE)
# lage tabell via excel:
#write_xlsx(felles_slekter_fraksjon, "/Users/siljerosbak/Dropbox/Master - Silje/Data/felles_slekter_fraksjon.xlsx")




# se hvor ulike slekter finnes
## noen i de pot giftige:
marinoscillum <- filter(Genus_unique_2_df, Genus == "Marinoscillum")
kordia <- filter(Genus_unique_2_df, Genus == "Kordia")
## noen i de ikke giftige:
oceanicaulis <- filter(Genus_unique_2_df, Genus == "Oceanicaulis")
Marivirga <- filter(Genus_unique_2_df, Genus == "Marivirga")
oiwenweeksia <- filter(Genus_unique_2_df, Genus == "Owenweeksia")


