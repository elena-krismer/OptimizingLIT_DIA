library(tidyverse)
library(ggtext)

# --------------- Figure 1  ----------------------------------------------------
# ---- FDR confidence estimation using a LIT for low-input experiments. --------
# The number of detected peptides when searching LIT-DIA datasets from 1 ng and
# 10 ng protein on column with different scan settings with windows that are
# matched to maintain approximately a 2-second cycle time. Searches were
# performed using the entrapment approach, where in addition to searching a H.
# sapiens database, peptides from E. coli, C. elegans, and S. cerevisiae were
# also considered. The number of misidentified entrapment peptides, as well as
# the ratio of entrapment to total detected peptides are shown above bars.

path  <- "data/figure1/"
theme_set(theme_bw())

# load data --------------------------------------------------------------------

#Turbo
HECS_turbo_1 <- read.delim2(paste0(path,
  "20221213_084930_HECS_SN17_directDIA_LITDIA_MS1_Zooom_MS2_Turbo_HeLa_01ng_Report.tsv"))
HECS_turbo_10 <- read.delim2(paste0(path,
  "20221213_085906_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Turbo_HeLa_10ng_Report.tsv"))

#Rapid
HECS_rapid_1 <- read.delim2(paste0(path,
  "20221213_085219_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Rapid_HeLa_01ng_Report.tsv"))
HECS_rapid_10 <- read.delim2(paste0(path,
  "20221213_084452_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Rapid_HeLa_10ng_Report.tsv"))

#Normal
HECS_normal_1 <- read.delim2(paste0(path,
  "20221213_083607_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Normal_HeLa_01ng_Report.tsv"))
HECS_normal_10 <- read.delim2(paste0(path,
  "20221213_084640_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Normal_HeLa_10ng_Report.tsv"))

#Enhanced
HECS_enhanced_1 <- read.delim2(paste0(path,
  "20221213_085737_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Enhanced_HeLa_01ng_Report.tsv"))
HECS_enhanced_10 <- read.delim2(paste0(path,
  "20221213_084822_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Enhanced_HeLa_10ng_Report.tsv"))

#Zoom
HECS_zoom_1 <- read.delim2(paste0(path,
 "20221213_090214_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Zoom_HeLa_01ng_Report.tsv"))
HECS_zoom_10 <- read.delim2(paste0(path,
"20221213_090106_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Zoom_HeLa_10ng_Report.tsv"))



HECS_turbo_1$input <- HECS_rapid_1$input <- HECS_normal_1$input <- HECS_enhanced_1$input <- HECS_zoom_1$input <- 1
HECS_turbo_10$input <- HECS_rapid_10$input <- HECS_normal_10$input <- HECS_enhanced_10$input <- HECS_zoom_10$input <- 10

HECS_turbo_1$MS2_speed <- HECS_turbo_10$MS2_speed <- "MS2-Turbo, 80W"
HECS_rapid_1$MS2_speed <- HECS_rapid_10$MS2_speed <- "MS2-Rapid, 57W"
HECS_normal_1$MS2_speed <- HECS_normal_10$MS2_speed <- "MS2-Normal, 40W"
HECS_enhanced_1$MS2_speed <- HECS_enhanced_10$MS2_speed <- "MS2-Enhanced, 25W"
HECS_zoom_1$MS2_speed <- HECS_zoom_10$MS2_speed <- "MS2-Zoom, 5W"

columns_to_extract <- c("PEP.StrippedSequence", "R.Replicate", "input", "MS2_speed", "PG.FastaFiles")

prot <-
  rbind(HECS_turbo_1[,columns_to_extract ],
        HECS_turbo_10[,columns_to_extract],
        HECS_rapid_1[, columns_to_extract],
        HECS_rapid_10[, columns_to_extract],
        HECS_normal_1[, columns_to_extract],
        HECS_normal_10[, columns_to_extract],
        HECS_enhanced_1[, columns_to_extract],
        HECS_enhanced_10[, columns_to_extract],
        HECS_zoom_1[, columns_to_extract],
        HECS_zoom_10[, columns_to_extract])


# clean data --------------------------------------------------------------------
# proteins  attributed to several organisms are removed.
# Names of the organisms are cleaned.

prot <- prot[!grepl(";", prot$PG.FastaFiles),]
prot$PG.FastaFiles %>% table()

prot$PG.FastaFiles <- gsub("uniprot-homo\\+sapiens-filtered-reviewed_yes", "Homo sapiens", prot$PG.FastaFiles)
prot$PG.FastaFiles <- gsub("uniprot-c\\+elegans-filtered-reviewed_yes", "C. elegans", prot$PG.FastaFiles)
prot$PG.FastaFiles <- gsub("uniprot-ecoli", "E. coli", prot$PG.FastaFiles)
prot$PG.FastaFiles <- gsub("uniprot_yeast", "Yeast", prot$PG.FastaFiles)

prot$PG.FastaFiles %>% table()

prot$PG.FastaFiles <- factor(prot$PG.FastaFiles, levels = c("E. coli", "Yeast", "C. elegans", "Homo sapiens"))
prot$MS2_speed <- factor(prot$MS2_speed, levels = unique(prot$MS2_speed))




# compute FDR ------------------------------------------------------------------

prot$positive <- prot$PG.FastaFiles

prot$positive <- sub("Homo sapiens", TRUE, prot$positive)
prot$positive <- sub("E. coli|C. elegans|Yeast", FALSE, prot$positive)


# The number of true positives and false positives peptides is counted.
# FDR is calculated by dividing false positives number by total positives (all peptides)


prot_tp_fp_count <- prot %>%
  group_by(R.Replicate, input, MS2_speed, positive) %>%
  summarize(pep_number = length(unique(PEP.StrippedSequence))) %>%
  group_by(input, MS2_speed, positive) %>%
  summarise(mean_pep = mean(pep_number, na.rm = TRUE))


fdr <-
  round(
    (1.152*(prot_tp_fp_count %>%
              filter(positive == FALSE) %>%
              pull(mean_pep))/
       prot_tp_fp_count %>% filter(positive == TRUE) %>%
       pull(mean_pep)) *100, 2)



prot_count$y_label <- prot_count$mean_pep

prot_count$y_label[prot_count$PG.FastaFiles == "E. coli"] <-
  prot_count$y_label[prot_count$PG.FastaFiles  == "E. coli"] +
  prot_count$y_label[prot_count$PG.FastaFiles  == "Yeast"] +
  prot_count$y_label[prot_count$PG.FastaFiles  == "C. elegans"] +
  prot_count$y_label[prot_count$PG.FastaFiles  == "Homo sapiens"] + 790

prot_count$y_label[prot_count$PG.FastaFiles == "Yeast"] <-
  prot_count$y_label[prot_count$PG.FastaFiles  == "Yeast"] +
  prot_count$y_label[prot_count$PG.FastaFiles  == "C. elegans"] +
  prot_count$y_label[prot_count$PG.FastaFiles == "Homo sapiens"] + 520

prot_count$y_label[prot_count$PG.FastaFiles  == "C. elegans"] <-
  prot_count$y_label[prot_count$PG.FastaFiles  == "C. elegans"] +
  prot_count$y_label[prot_count$PG.FastaFiles  == "Homo sapiens"] + 250

prot_count$y_label[prot_count$PG.FastaFiles  == "Homo sapiens"] <- NA

# Plot -------------------------------------------------------------------------

new_label <- c(rep("misidentified\nn = 10", 3), "Homo sapiens",rep("misidentified\nn = 9", 3), "Homo sapiens",
               rep("misidentified\nn = 14", 3), "Homo sapiens",rep("misidentified\nn = 14", 3), "Homo sapiens",
               rep("misidentified\nn = 2", 3), "Homo sapiens", rep("misidentified\nn = 24", 3), "Homo sapiens",
               rep("misidentified\nn = 13", 3), "Homo sapiens", rep("misidentified\nn = 17", 3), "Homo sapiens",
               rep("misidentified\nn = 8", 3), "Homo sapiens", rep("misidentified\nn = 5", 3), "Homo sapiens"
)
new_y_pos <- c(rep(2794, 3), NA,
                 rep(3028, 3), NA,
                 rep(3345, 3), NA,
                 rep(2558, 3), NA,
                 rep(1089, 3), NA,
                 rep(5909,3), NA,
                 rep(5423, 3), NA,
                 rep(5086, 3), NA,
                 rep(3311, 3), NA,
                 rep(1287, 3), NA) + 100


fdr_label <- c(rep(fdr[1], 4), rep(fdr[2], 4), rep(fdr[3], 4), rep(fdr[4], 4), rep(fdr[5], 4),
                 rep(fdr[6], 4), rep(fdr[7], 4), rep(fdr[8], 4), rep(fdr[9], 4), rep(fdr[10], 4))
fdr_label <- paste0(fdr_label_2, "%")
fdr_label_pos<- rep(c(6900, NA, NA, NA), 10)

figure_1 <- prot_count %>%
  ggplot(aes(x = as.factor(input), y = mean_pep,
             fill = organism)) +
  geom_col(color = "black", size = 0.08, alpha=0.8)+
  facet_grid(~MS2_speed)+
  geom_label(aes(label = new_label, y = new_y_pos),
             size = 4, show.legend = FALSE, alpha=0.8) +
  scale_fill_manual(values = c( "#009599", "grey87"))+
  labs(x = "Input (ng)",
       y = "Identified peptides",
       fill = "Organism")+
  scale_y_continuous(limits = c(0, 6900))+
  theme(axis.title = element_text(size = axis_title),
        axis.text = element_text(size = axis_text),
        strip.text = element_text(size=16),
        strip.background = element_rect(fill = "#E5E4F7"),
        legend.text = element_markdown(size=axis_text),
        legend.title = element_text(size= axis_title),
        legend.position = "top")+
  geom_text(aes(label = fdr_label, y = fdr_label_pos), size = 4, show.legend = FALSE)



