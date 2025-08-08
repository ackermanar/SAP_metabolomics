library(data.table)
library(pheatmap)
library(grid)
library(gridExtra)
library(ppcor)
library(tidyverse)
library(wesanderson)
library(RColorBrewer)
library(viridis)
library(gridExtra)
library(ggforce)
library(EnhancedVolcano)

setwd("SAP_Metabolomics")

# Figure 1 - RACE

met <- fread("./data/raw/untreated_CU_met/SAP19_MET.csv")
data1 <- read.csv("./data/processed/Fig1/Data/final_peak_list_RACE.csv") #with controls

filtered1 <- data1 %>%
  filter(Hits.sig > 2) %>%
  mutate(RACE = ifelse(RACE == "Tx2911", "Tx2911 (Check Line)", RACE)) %>%
  mutate(RACE = ifelse(RACE == "P850029", "P850029 (Check Line)", RACE))

race_counts <- met %>%
  group_by(RACE) %>%
  summarize(Count = n())%>%
  mutate(RACE = ifelse(RACE == "Tx2911", "Tx2911 (Check Line)", RACE)) %>%
  mutate(RACE = ifelse(RACE == "P850029", "P850029 (Check Line)", RACE))

race_counts2 <- filtered1 %>%
  group_by(RACE) %>%
  summarize()%>%
  mutate(RACE = ifelse(RACE == "Tx2911", "Tx2911 (Check Line)", RACE)) %>%
  mutate(RACE = ifelse(RACE == "P850029", "P850029 (Check Line)", RACE))

merged_counts1 <- merge(race_counts2, race_counts, by = "RACE", all.x = TRUE)
merged_counts1$Label <- paste(merged_counts1$RACE, "\n (n = ", merged_counts1$Count, ")", sep = "")

ggplot(filtered1, aes(x = RACE, y = value, size = EASE_log10, color = EASE_log10)) +
  geom_point() +
  ylab("Pathway") +
  xlab("Race") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = merged_counts1$Label) +
  scale_size_continuous(limits = c(0, 3), breaks = c(0, 1, 2, 3), range = c(1, 10)) +
  scale_color_gradient2(limits = c(0, 3), breaks = c(0, 1, 2, 3), low = "lightblue1", mid = "lightskyblue", high = "dodgerblue4") +
  guides(size = guide_legend(), color = guide_legend()) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 16))
ggsave("Fig1_Race.tiff", device = "tiff", width = 14, height = 12, dpi = 300)

# Figure 2 - PANICLE STRUCTURE
met <- fread("./data/processed/Fig2/Data/final_peak_list_PANICLE.csv")
data2 <- read.csv("./data/processed/Fig2/Data/final_peak_list_PANICLE.csv")

filtered2 <- data2 %>%
  filter(Hits.sig > 2)

panicle_counts <- met %>%
  group_by(PANICLE) %>%
  summarize(Count = n())

panicle_counts2 <- filtered2 %>%
  group_by(PANICLE) %>%
  summarize()

merged_counts2 <- merge(panicle_counts2, panicle_counts, by = "PANICLE", all.x = TRUE)
merged_counts2$Label <- paste(merged_counts2$PANICLE, "\n (n = ", merged_counts2$Count, ")", sep = "")

ggplot(filtered2, aes(x = PANICLE, y = value, size = EASE_log10, color = EASE_log10)) +
  geom_point() +
  ylab("Pathway") +
  xlab("Panicle Structure") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = merged_counts2$Label) +
  scale_size_continuous(limits = c(0, 3), breaks = c(0, 1, 2, 3), range = c(5, 15)) +
  scale_color_gradient2(limits = c(0, 3), breaks = c(0, 1, 2, 3), low = "lightblue1", mid = "lightskyblue", high = "dodgerblue4") +
  guides(color = guide_legend(), size = guide_legend()) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16))

#Large Print
ggsave("Fig2_Pan.tiff", device = "tiff", width = 14, height = 12, dpi = 300)

# Figure 3 - GRAIN COLOR

data3 <- read_csv("./data/processed/Fig3/Data/final_peak_list_COLOR.csv")

filtered3 <- data3 %>%
  filter(Hits.sig>2)

pericarp_counts <- met %>% 
  group_by(PERICARP) %>%
  summarize(Count = n())

pericarp_counts2 <- filtered3 %>% 
  group_by(PERICARP) %>%
  summarize()

merged_counts3 <- merge(pericarp_counts2, pericarp_counts, by = "PERICARP", all.x = TRUE)
merged_counts3$Label <- paste(merged_counts3$PERICARP, "\n (n = ", merged_counts3$Count, ")", sep = "")

ggplot(filtered3, aes(x = PERICARP, y = value, size = EASE_log10, color = EASE_log10)) +
  geom_point() +
  ylab("Pathway") +
  xlab("Grain Color") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = merged_counts3$Label) +
  scale_size_continuous(limits = c(0, 3), breaks = c(0, 1, 2, 3), range = c(3, 12)) +
  scale_color_gradient2(limits = c(0, 3), breaks = c(0, 1, 2, 3), low = "lightblue1", mid = "lightskyblue", high = "dodgerblue4") +
  guides(color = guide_legend(), size = guide_legend()) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16))


#Large Print
ggsave("Fig3_Per.tiff", device = "tiff", width = 12, height = 14, dpi = 300)

# Figure 4
sgmScores <- read_csv("./data/raw/field_data/SGM_DATA.csv")

ggplot(sgmScores, aes(x = PGMSR, y = Fumonisin, , size = FSDI, color = Geno)) +
  geom_point(alpha = .75) +
  ylab("Fumonisin (ppm)") +
  xlab("PGMSR") +
  scale_color_manual(values = c("Tx2911" = "darkorange2", "P850029" = "dodgerblue2"), name = "Accession") +  scale_size(range = c(4, 20)) +
  theme(legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        axis.title = element_text(size = 16))

#Large Print
ggsave("Fig4_SMG_large.tiff", device = "tiff", width = 14, height = 12, dpi = 300)
#Smaller Print
ggsave("Fig4_SGM_small.tiff", device = "tiff", width = 7, height = 6, dpi = 300)

# Figure 5
pcaScores <- read_csv("./data/processed/Fig5/Data/pca_score.csv")

ggplot(pcaScores, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(aes(fill = after_scale(alpha(color, 0.75))), shape = 21, size = 12) +
  stat_ellipse(geom = 'polygon', alpha = .1, aes(fill = Group)) +
  labs(x = "PC1", y = "PC2") +
  scale_color_manual(values = c("dodgerblue2", "darkorange2"), name = "Accession") +
  scale_fill_manual(values = c("dodgerblue4", "darkorange4")) +
  theme(legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  guides(fill = "none", alpha = "none", size = "none")

#Large Print
ggsave("Fig5_PCA_large.tiff", device = "tiff", width = 14, height = 12, dpi = 300)
#Smaller Print
ggsave("Fig5_PCA_small.tiff", device = "tiff", width = 7, height = 6, dpi = 300)

# Figure 6
ftAnn <- read_csv("./data/raw/annotatesummary.csv")
lfc2 <- read_csv("./data/processed/Fig7/Data/Tx2911vol.csv")
full <- read_csv("./data/processed/Fig7/Data/Tx2911upreg.csv")

sumFt <- ftAnn %>%
  left_join(lfc2, by = "alignment_id") %>%
  filter(!is.na(log2FC)) %>%
  group_by(Class) %>%
  summarise(count = n(), log2FoldChange = mean(log2FC), SE = sd(log2FC)/sqrt(n()), SD = sd(log2FC))

sumFt <- full %>%
    filter(!is.na(log2FoldChange)) %>%
    group_by(Class) %>%
    filter(Class != "Unknown") %>%
    summarise(count = n(), log2FoldChangeAVG = mean(log2FoldChange), SE = sd(log2FoldChange)/sqrt(n()), SD = sd(log2FoldChange)) %>%
    mutate(label = paste0("n = ", count)) %>%
    mutate(label = paste0("italic('", label, "')")) %>%
    arrange(log2FoldChangeAVG)

ggplot(sumFt, aes(x=reorder(Class, log2FoldChangeAVG), y=log2FoldChangeAVG, fill = log2FoldChangeAVG)) +
  geom_bar(width = 1, stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChangeAVG - SE, ymax = log2FoldChangeAVG + SE), width = 0.2, color = "cornsilk4") +
  coord_flip() +
  scale_fill_gradient2(limits = c(-2.35, 4), breaks = c(-2.35, 0, 2, 4), low = "dodgerblue2", mid = "cyan3", high = "darkorange2") +
  geom_text(aes(label = label), parse = TRUE, vjust = -0.5, size = 5) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))  +
  theme(legend.position = "none") +
  labs(title = "Overview of Annotated Features in P850029 and Tx2911", size = 16) +
  ylab(expression(paste("P850029                              Average ", Log[2], " Fold Change                           Tx2911"))) +
  theme(legend.text = element_text(size = 12),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.title.x = element_text(size = 16, hjust = 0),
        axis.title.y=element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5))

ggsave("Fig6_2SidedBar.tiff", width = 12, height = 16, dpi = 300, units = "in")

# Figure 7

res <- read.csv("./data/processed/Fig7/Data/Tx2911upreg.csv", row.names = 1)

keyvals.colour <- ifelse(
    res$log2FoldChange < -1.25, "dodgerblue2",
      ifelse(res$log2FoldChange > 1.25, "darkorange2",
        "ivory3"))
  keyvals.colour[is.na(keyvals.colour)] <- "black"
  names(keyvals.colour)[keyvals.colour == "darkorange2"] <- 'Abundant in Tx2911'
  names(keyvals.colour)[keyvals.colour == "ivory4"] <- 'Shared by both'
  names(keyvals.colour)[keyvals.colour == "dodgerblue2"] <- 'Abundant in P850029'
  
  keyvals.shape <- ifelse(
    res$Class == 'Flavonoids', 19,
      ifelse(res$Class == 'Phenylpropanoids', 17,
                    ifelse(res$Class == 'Phenolic_acids', 18,
                           ifelse(res$Class == 'Phenols', 15,
                                  ifelse(res$Class == 'Glycerophospholipids', 4,
        1)))))
  keyvals.shape[is.na(keyvals.shape)] <- 1
 
  names(keyvals.shape)[keyvals.shape == 19] <- 'Flavonoids'
  names(keyvals.shape)[keyvals.shape == 17] <- 'Phenylpropanoids'
  names(keyvals.shape)[keyvals.shape == 18] <- 'Phenolic acids'
  names(keyvals.shape)[keyvals.shape == 15] <- 'Phenols'
  names(keyvals.shape)[keyvals.shape == 4] <- 'Glycerophospholipids'
  names(keyvals.shape)[keyvals.shape == 1] <- 'Other/Unknown'



vol <- EnhancedVolcano(res,
                       lab = rownames(res),
                       x = 'log2FoldChange',
                       y = 'pvalue',
                       xlim = c(-6, 11),
                       title = "P850029 vs Tx2911 ",
                       pCutoff = .05,
                       selectLab = c('Luteolin','Apigeninidin', 'Caffeic_Acid', 'Salicylic_acid', 'Luteolin', 'Luteolinidin', 'Naringenin', 'Glycerylphosphorylcholine'),
                       drawConnectors = TRUE,
                       FCcutoff = 1.25,
                       pointSize = 5.0,
                       shapeCustom = keyvals.shape,
                       colCustom = keyvals.colour,
                       colAlpha = 5/5
                       )

print(vol)

ggsave("./data/processed/Fig7/Fig7_Volcano.tiff", plot = vol, width = 16, height = 12, dpi = 300, units = "in")

# Figure 8: Heatmap of metabolites and samples
data <- read_csv("./data/processed/Fig8/Data/data_normalized.csv") %>%
  column_to_rownames("Sample")
data <- t(data)

metadata <- read.csv("./data/processed/Fig8/Data/metadata.csv", row.names = 1)

rowData <- read_csv("./data/processed/Fig8/Data/correlation_feature.csv") %>%
    column_to_rownames("feature") %>%
    dplyr::select(sig, correlation) %>%
    rename(Correlation = correlation, P.value = sig)

rowOrder <- c("Glutathione", "Quercetin 4'-O-beta-D-glucopyranoside", "Corymboside", "Kynurenic Acid", "PubChem CID 6426860*", "Vanillic Acid", "Caffeic Acid", "Chlorogenic Acid", "Luteolin", "Luteolinidin", "Naringenin", "Salicylic Acid", "Vanillin", "Apigenidin")

colOrder <- c("P85_Rep1-188", "P85_Rep1-210", "P85_Rep1-258", "P85_Rep1-429", "P85_Rep1-440", "P85_Rep1-521", "P85_Rep1-550", "P85_Rep1-552", "P85_Rep1-602", "P85_Rep1-613", "P85_Rep1-639", "P85_Rep2-11", "P85_Rep2-123", "P85_Rep2-134", "P85_Rep2-151", "P85_Rep2-378", "P85_Rep2-547", "P85_Rep2-549", "P85_Rep2-565", "P85_Rep2-68", "Tx2911_Rep1-308", "Tx2911_Rep1-333", "Tx2911_Rep1-468", "Tx2911_Rep1-542", "Tx2911_Rep1-543", "Tx2911_Rep1-676", "Tx2911_Rep2-101", "Tx2911_Rep2-181", "Tx2911_Rep2-205", "Tx2911_Rep2-345", "Tx2911_Rep2-392", "Tx2911_Rep2-512", "Tx2911_Rep2-529")

colors <- brewer.pal(9, "Blues")
ann_colors <- list(
    Label = c("P85" = "dodgerblue2", "Tx2911" = "darkorange2"),
    P.value = c("Significant" = "chartreuse2", "Insignificant" = "coral2"),
    Correlation = viridis(100)
)
filename = "hmp2.tiff"
tiff(filename, width = 16, height = 12, units = "in", res = 300, type = "cairo")
pheatmap(data[rowOrder, colOrder],
    annotation_col = metadata,
    annotation_row = rowData,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    cellheight = 30,
    cellwidth = 15,
    fontsize_row = 12,
    gaps_col = 20,
    width = 16,
    annotation_colors = ann_colors,
    annotation_names_row = FALSE,
    annotation_names_col = TRUE,
    show_colnames = FALSE,
    show_rownames = FALSE,
    filename = filename)

# Supplementary Figure 1 - COLOR and PAN.STRUCTURE across RACE

met2 <- met %>%
  group_by(PERICARP, RACE) %>%
  na.omit() %>%
  summarize(count = n())

merged_pericarp_race <- merge(pericarp_counts2, met2, by = "PERICARP", all.x = TRUE)
merged_pericarp_race2 <- merge(race_counts2, merged_pericarp_race, by = "RACE", all.x = FALSE)

met3 <- met %>%
  group_by(PANICLE, RACE) %>%
  na.omit() %>%
  summarize(count = n())

merged_panicle_race <- merge(panicle_counts2, met3, by = "PANICLE", all.x = TRUE)
merged_panicle_race2 <- merge(race_counts2, merged_panicle_race, by = "RACE", all.x = FALSE)

race_per <- ggplot(merged_pericarp_race2, aes(x = RACE, y = PERICARP, size = count, color = count)) +
  geom_point() +
  ylab("Grain Color") +
  xlab("Race") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())+
  scale_size_continuous(limits = c(1, 100), breaks = c(1, 25, 50, 75, 100), range = c(0, 12)) +
  scale_color_gradient2(limits = c(1, 100), breaks = c(1, 25, 50, 75, 100), low = "lightblue1", mid = "lightskyblue", high = "dodgerblue4") +
  guides(color = guide_legend(), size = guide_legend()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 16)) +
  theme(legend.position = "top")

race_pan <- ggplot(merged_panicle_race2, aes(x = RACE, y = PANICLE, size = count, color = count)) +
  geom_point() +
  ylab("Panicle Structure") +
  xlab("Race") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size_continuous(limits = c(1, 100), breaks = c(1, 25, 50, 75, 100), range = c(0, 12)) +
  scale_color_gradient2(limits = c(1, 100), breaks = c(1, 25, 50, 75, 100), low = "lightblue1", mid = "lightskyblue", high = "dodgerblue4") +
  guides(size = guide_legend(title = "Count"), color = guide_legend(title = "Count")) +
  theme(legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  theme(legend.position = "none")

combined_plot <- ggarrange(race_per, race_pan, ncol = 1, align = "v")
tiff("combined_plot.tiff", width = 10, height = 10, units = "in", res = 300)
grid.draw(combined_plot)
dev.off()
