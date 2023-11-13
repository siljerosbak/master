# Statistikk

# Pakker jeg trenger
library(phyloseq)
library(ggplot2)   
library(ggpubr)
library(readxl)       # necessary to import the data from Excel file
library(dplyr)        # filter and reformat data frames
library(tibble)       # Needed for converting column to row names
library(knitr)
library(tidyverse)
library(vegan)
library(cowplot)


# Laster inn filene mine og setter wd:

knitr::opts_knit$set(root.dir = "/Users/siljerosbak/Dropbox/Master - Silje/Data/") # hvor jeg henter dataene fra

metadata <- read_excel("Metadata1_master.xlsx") # laster inn metadataene
samples_df <- readRDS("samples_df.rds", refhook=NULL) # metadata
rarefied_ps <- readRDS("rarefied_ps.rds", refhook=NULL) # rarefied phyloseq 
p_rarefied_df <- readRDS("p_rarefied_df.rds", refhook=NULL)
p_clean_gg_uten_2 <- readRDS("p_clean_gg_uten_2.rds", refhook=NULL) # phyloseq uten de negative kontroll filterene mine


# legger til farger
color_palette_class_2 = c("#D1BBD7", "#AE76A3", "#882E72", "#1965B0",
                          "#5289C7", "#7BAFDE", "#4EB265", "#90C987",
                          "#CAE0AB", "#F7F056", "#F6C141", "#F1932D",
                          "#E8601C", "#DC050C", "#72190E")



# Tabell med gjennomsnittlig alfadiversitet for alle prøvene:

alpha_diversity <- estimate_richness(rarefied_ps, measures = c("Shannon", "Simpson", "Chao1", "ACE", "Fisher"))
alpha_div_matrix <- as.data.frame(alpha_diversity)

alpha_diversity_df <- data.frame("Prøve"=rownames(alpha_div_matrix), alpha_div_matrix)

plot_richness(rarefied_ps, x="Label", measures=c("Shannon", "Simpson", "ACE"), color="Giftighet") +
  labs(x = "Prøvenummer",
       y = bquote(alpha~"Diversitet")) + 
  scale_color_manual(values=c("#EF7F4FFF", "#90C987")) +
  geom_point(aes(color=Giftighet)) +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=8),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        panel.grid.major = element_line(color = "lightgray"),
        panel.grid.minor = element_line(color = "lightgray"),
        panel.border = element_rect(color = "gray", fill = NA, size = .5)) 



# T-test: for å sammenligne ikke-giftig og giftig
# Ordner på dataene slik at jeg kan merge:
replace(rownames(samples_df), "-", ".")
ny_alpha_div_matrix <- alpha_div_matrix[order(rownames(alpha_div_matrix)), ]
ny_samples_df <- samples_df[order(rownames(samples_df)), ]
samples_df_alpha<- ny_samples_df[-c(4,21),] # Fjerner negative

# Merge the data frames based on SampleID
merged_alpha <- cbind(ny_alpha_div_matrix, samples_df_alpha)

merged_alpha$Giftighet <- as.factor(merged_alpha$Giftighet)

# Sjekker assumtpions 
par(mfrow = c(2, 2))
ggplot(merged_alpha, aes(x = ACE)) +
  geom_histogram(aes(color = Giftighet, fill = Giftighet), 
                 position = "identity", bins = 30, alpha = 0.4) +
  scale_color_manual(values = c("#F7F056", "#F1932D")) +
  scale_fill_manual(values = c("#F7F056", "#F1932D"))

ggplot(merged_alpha, aes(x = Simpson)) +
  geom_histogram(aes(color = Giftighet, fill = Giftighet), 
                 position = "identity", bins = 30, alpha = 0.4) +
  scale_color_manual(values = c("#F7F056", "#F1932D")) +
  scale_fill_manual(values = c("#F7F056", "#F1932D"))

ggplot(merged_alpha, aes(x = Shannon)) +
  geom_histogram(aes(color = Giftighet, fill = Giftighet), 
                 position = "identity", bins = 30, alpha = 0.4) +
  scale_color_manual(values = c("#F7F056", "#F1932D")) +
  scale_fill_manual(values = c("#F7F056", "#F1932D"))

ggplot(merged_alpha, aes(x = Fisher)) +
  geom_histogram(aes(color = Giftighet, fill = Giftighet), 
                 position = "identity", bins = 30, alpha = 0.4) +
  scale_color_manual(values = c("#F7F056", "#F1932D")) +
  scale_fill_manual(values = c("#F7F056", "#F1932D"))


#gjør welchs test på de ikke-normalfordelte + liten sample size
t.test(Simpson ~ Giftighet, data = merged_alpha)
t.test(Shannon ~ Giftighet, data = merged_alpha)
t.test(ACE ~ Giftighet, data = merged_alpha)
t.test(Fisher ~ Giftighet, data = merged_alpha)
t.test(Chao1 ~ Giftighet, data = merged_alpha) 



# Box-plot med alfa-diversiteter
a1 <- ggboxplot(merged_alpha, x = "Giftighet", y = "ACE",
                add = "dotplot",
                color = "Giftighet", palette =c("#F7F056", "#F1932D"),) +
  (stat_compare_means(method = "t.test")) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) + #remove x axis labels
  theme(text = element_text(size = 18))


a2 <- ggboxplot(merged_alpha, x = "Giftighet", y = "Shannon",
                add = "dotplot",
                color = "Giftighet", palette =c("#F7F056", "#F1932D"),) +
  (stat_compare_means(method = "t.test")) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) + #remove x axis labels
  theme(text = element_text(size = 18))



a3 <- ggboxplot(merged_alpha, x = "Giftighet", y = "Simpson",
                add = "dotplot",
                color = "Giftighet", palette =c("#F7F056", "#F1932D"),) +
  (stat_compare_means(method = "t.test")) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) + #remove x axis labels
  theme(text = element_text(size = 18))


a4 <- ggboxplot(merged_alpha, x = "Giftighet", y = "Fisher",
                add = "dotplot",
                color = "Giftighet", palette =c("#F7F056", "#F1932D"),) +
  (stat_compare_means(method = "t.test")) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) + #remove x axis labels
  theme(text = element_text(size = 18))

plot_grid(a1, a2, a3, a4,
          labels = c('A', 'B', 'C', 'D'),
          align="none")  






#Plot som kan begrunner hvorfor de to fraksjonene er forskjellige
merged_alpha$id <- with(merged_alpha, paste(Stamme...4, Sampling)) # legger til en ID som forteller stemma + hvilket "replikat"

attached <- filter(merged_alpha, Filter_strl_um == "Fastsittende(1um)")

freeliving <- filter(merged_alpha, Filter_strl_um == "Frittlevende(0,2um)")

s <- merge(attached, freeliving, by= "id") #slår de sammen igjen slik at alt er i en df

ggplot(s, aes(x=Shannon.x, y=Shannon.y, color=Giftighet.x)) + 
  geom_point(size=6)

ggplot(s, aes(x=Chao1.x, y=Chao1.y, color=Giftighet.x)) + 
  geom_point(size=6) +
  geom_abline(intercept=0, slope=1)
#theme_ipsum()

#Lager et plot som viser relasjonen mellom alle alfa-diversitetsmålene 
plot(exp(Shannon) ~ ACE, data=merged_alpha)
plot((Shannon) ~ ACE, data=merged_alpha)
plot(merged_alpha[, c(1,3,5,6,7)]) # lager et felles plot for de ulike alfadiversitets-målene som viser sammenhengen mellom dem

# Merge the data frames based on SampleID
merged_alpha <- cbind(ny_alpha_div_matrix, samples_df_alpha)

merged_alpha$Filter_strl_um <- as.factor(merged_alpha$Filter_strl_um)

# Sjekker assumptions
par(mfrow = c(2, 2))
ggplot(merged_alpha, aes(x = ACE)) +
  geom_histogram(aes(color = Filter_strl_um, fill = Filter_strl_um), 
                 position = "identity", bins = 30, alpha = 0.4) +
  scale_color_manual(values = c("#F7F056", "#F1932D")) +
  scale_fill_manual(values = c("#F7F056", "#F1932D"))

ggplot(merged_alpha, aes(x = Simpson)) +
  geom_histogram(aes(color = Filter_strl_um, fill = Filter_strl_um), 
                 position = "identity", bins = 30, alpha = 0.4) +
  scale_color_manual(values = c("#F7F056", "#F1932D")) +
  scale_fill_manual(values = c("#F7F056", "#F1932D"))

ggplot(merged_alpha, aes(x = Shannon)) +
  geom_histogram(aes(color = Filter_strl_um, fill = Filter_strl_um), 
                 position = "identity", bins = 30, alpha = 0.4) +
  scale_color_manual(values = c("#F7F056", "#F1932D")) +
  scale_fill_manual(values = c("#F7F056", "#F1932D"))

ggplot(merged_alpha, aes(x = Fisher)) +
  geom_histogram(aes(color = Filter_strl_um, fill = Filter_strl_um), 
                 position = "identity", bins = 30, alpha = 0.4) +
  scale_color_manual(values = c("#F7F056", "#F1932D")) +
  scale_fill_manual(values = c("#F7F056", "#F1932D"))

# Fordelingen ser lik ut (den er ikke normalfordelt, men den har antydning til samme fordeling)




# t-test for filterstørls
t.test(Simpson ~ Filter_strl_um, data = merged_alpha)
t.test(Shannon ~ Filter_strl_um, data = merged_alpha)
t.test(ACE ~ Filter_strl_um, data = merged_alpha)
t.test(Fisher ~ Filter_strl_um, data = merged_alpha) 
t.test(Chao1 ~ Filter_strl_um, data = merged_alpha) 

# Box-plot med alfa-diversiteter
a1 <- ggboxplot(merged_alpha, x = "Filter_strl_um", y = "ACE",
                add = "dotplot",
                color = "Filter_strl_um", palette =c("#F7F056", "#F1932D"),) +
  (stat_compare_means(method = "t.test")) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) + #remove x axis labels
  theme(text = element_text(size = 18))


a2 <- ggboxplot(merged_alpha, x = "Filter_strl_um", y = "Shannon",
                add = "dotplot",
                color = "Filter_strl_um", palette =c("#F7F056", "#F1932D"),) +
  (stat_compare_means(method = "t.test")) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) + #remove x axis labels
  theme(text = element_text(size = 18))



a3 <- ggboxplot(merged_alpha, x = "Filter_strl_um", y = "Simpson",
                add = "dotplot",
                color = "Filter_strl_um", palette =c("#F7F056", "#F1932D"),) +
  (stat_compare_means(method = "t.test")) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) + #remove x axis labels
  theme(text = element_text(size = 18))


a4 <- ggboxplot(merged_alpha, x = "Filter_strl_um", y = "Fisher",
                add = "dotplot",
                color = "Filter_strl_um", palette =c("#F7F056", "#F1932D"),) +
  (stat_compare_means(method = "t.test")) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank()) + #remove x axis labels
  theme(text = element_text(size = 18))

plot_grid(a1, a2, a3, a4,
          labels = c('A', 'B', 'C', 'D'),
          align="none")       




# Beta diversity 
## NMDS plot
ps.prop <- transform_sample_counts(p_clean_gg_uten_2, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="MDS", distance="bray")

NMDS_plot_1 <- plot_ordination(ps.prop, ord.nmds.bray, color="Art_stamme", shape = "Filter_strl_um") + 
  labs(x = "NMDS1",
       y = "NMDS2") +
  scale_color_manual(values=c("#AE76A3","#7BAFDE", "#90C987","#1965B0", 
                              "#F7F056", "#F1932D", "#E8601C", "#DC050C", "#72190E")) +
  geom_point(aes(color=Art_stamme), alpha = 0.7, size = 4) +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=8),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        panel.border = element_rect(color = "gray", fill = NA, size = .5)) +
  coord_fixed()

NMDS_plot_1
#ggsave("NMDS_1.pdf", dpi=600)

NMDS_plot_2 <- plot_ordination(ps.prop, ord.nmds.bray, color="Art_stamme", shape = "Giftighet") + 
  labs(x = "NMDS1",
       y = "NMDS2") +
  scale_color_manual(values=c("#AE76A3","#7BAFDE", "#90C987","#1965B0", 
                              "#F7F056", "#F1932D", "#E8601C", "#DC050C", "#72190E")) +
  geom_point(aes(color=Art_stamme), alpha = 0.7, size = 4) +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=8),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        panel.border = element_rect(color = "gray", fill = NA, size = .5)) +
  coord_fixed()

NMDS_plot_2
#ggsave("NMDS_2.pdf", dpi=600)


NMDS_plot_4 <- plot_ordination(ps.prop, ord.nmds.bray, color="Slekt", shape = "Giftighet") + 
  labs(x = "NMDS1",
       y = "NMDS2") +
  scale_color_manual(values=c( "#AE76A3",
                               "#7BAFDE", "#90C987",
                               "#CAE0AB", "#F7F056", "#F6C141", "#F1932D",
                               "#E8601C", "#DC050C", "#72190E")) +
  geom_point(aes(color=Slekt), alpha = 0.7, size = 4) +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=8),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        panel.border = element_rect(color = "gray", fill = NA, size = .5)) +
  coord_fixed()

NMDS_plot_4
#ggsave("NMDS_3.pdf", dpi=600)

NMDS_plot_5 <- plot_ordination(ps.prop, ord.nmds.bray, color="Slekt", shape = "Filter_strl_um") + 
  labs(x = "NMDS1",
       y = "NMDS2") +
  scale_color_manual(values=c( "#CAE0AB", "#F7F056", "#F1932D",
                               "#E8601C", "#DC050C", "#72190E")) +
  geom_point(aes(color=Slekt), alpha = 0.7, size = 4) +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=10, face="bold"),
        legend.text = element_text(size=8),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        panel.border = element_rect(color = "gray", fill = NA, size = .5)) +
  coord_fixed()

NMDS_plot_5
#ggsave("NMDS_4.png", dpi=600)


  


## PERMANOVA

dist.bray <- phyloseq::distance(rarefied_ps, method = "bray")
adonis_metadata <- as(sample_data(rarefied_ps), "data.frame")

#sjekker assumptions
distance1 <- vegdist(dist.bray, method = "bray")
dispersion1 <- betadisper(distance1, group = adonis_metadata$Giftighet)
permutest(dispersion1) #kan brukes om den er større enn 0.05 p-verdi (ikke-signifkant p-verdi)

#adonis/permanova
set.seed(100) 
adonis_1 = adonis2(dist.bray ~ Giftighet * Slekt, data = adonis_metadata, permutations=2000)
adonis_1

set.seed(100) 
adonis_2 = adonis2(dist.bray ~ Filter_strl_um * Slekt, data = adonis_metadata, permutations=2000)
adonis_2

