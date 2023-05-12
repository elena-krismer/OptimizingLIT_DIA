library(tidyverse)
library(ggtext)
library(readr)

# --------------- Figure 1  ----------------------------------------------------
# ---- FDR confidence estimation using a LIT for low-input experiments. --------
# The number of detected peptides when searching LIT-DIA datasets from 1 ng and
# 10 ng protein on column with different scan settings with windows that are
# matched to maintain approximately a 2-second cycle time. Searches were
# performed using the entrapment approach, where in addition to searching a H.
# sapiens database, peptides from E. coli, C. elegans, and S. cerevisiae were
# also considered. The number of misidentified entrapment peptides, as well as
# the ratio of entrapment to total detected peptides are shown above bars.

path  <- "data/figure1"
theme_set(theme_bw())

file_names <- c(
  "20221213_084930_HECS_SN17_directDIA_LITDIA_MS1_Zooom_MS2_Turbo_HeLa_01ng_Report.tsv",
  "20221213_085906_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Turbo_HeLa_10ng_Report.tsv",
  "20221213_085219_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Rapid_HeLa_01ng_Report.tsv",
  "20221213_084452_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Rapid_HeLa_10ng_Report.tsv",
  "20221213_083607_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Normal_HeLa_01ng_Report.tsv",
  "20221213_084640_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Normal_HeLa_10ng_Report.tsv",
  "20221213_085737_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Enhanced_HeLa_01ng_Report.tsv",
  "20221213_084822_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Enhanced_HeLa_10ng_Report.tsv",
  "20221213_090214_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Zoom_HeLa_01ng_Report.tsv",
  "20221213_090106_HECS_SN17_directDIA_LITDIA_MS1_Zoom_MS2_Zoom_HeLa_10ng_Report.tsv"
)

read_data <- function(file_path) {
  read_delim(file_path, delim = "\t", col_types = cols())
}

file_data <- file_names %>%
  set_names() %>%
  map_dfr(~read_data(file.path(path, .x))) %>%
  mutate(
    input = if_else(str_detect(file_names, "01ng"), 1, 10),
    MS2_speed = case_when(
      str_detect(file_names, "Turbo") ~ "MS2-Turbo, 80W",
      str_detect(file_names, "Rapid") ~ "MS2-Rapid, 57W",
      str_detect(file_names, "Normal") ~ "MS2-Normal, 40W",
      str_detect(file_names, "Enhanced") ~ "MS2-Enhanced, 25W",
      str_detect(file_names, "Zoom") ~ "MS2-Zoom, 5W"
    ),
    PG.FastaFiles = str_extract(file_names, "(?<=_HeLa_)\\w+")
  ) %>%
  select(PEP.StrippedSequence, R.Replicate, input, MS2_speed, PG.FastaFiles)

# combine data -----------------------------------------------------------------

columns_to_extract <- c("PEP.StrippedSequence", "R.Replicate", "input", "MS2_speed", "PG.FastaFiles")

prot <- file_data %>%
  filter(!str_detect(PG.FastaFiles, ";")) %>%
  select(all_of(columns_to_extract)) %>%
  mutate(
    PG.FastaFiles = str_replace(PG.FastaFiles, "uniprot-homo\\+sapiens-filtered-reviewed_yes", "Homo sapiens"),
    PG.FastaFiles = str_replace(PG.FastaFiles, "uniprot-c\\+elegans-filtered-reviewed_yes", "C. elegans"),
    PG.FastaFiles = str_replace(PG.FastaFiles, "uniprot-ecoli", "E. coli"),
    PG.FastaFiles = str_replace(PG.FastaFiles, "uniprot_yeast", "Yeast"),
    PG.FastaFiles = factor(PG.FastaFiles, levels = c("E. coli", "Yeast", "C. elegans", "Homo sapiens")),
    MS2_speed = factor(MS2_speed, levels = unique(MS2_speed))
  )


# compute FDR ------------------------------------------------------------------

prot$positive <- prot$PG.FastaFiles %>% str_detect("Homo sapiens")

prot$positive <- sub("Homo sapiens", TRUE, prot$positive)
prot$positive <- sub("E. coli|C. elegans|Yeast", FALSE, prot$positive)


# The number of true positives and false positives peptides is counted.
# FDR is calculated by dividing false positives number by total positives (all peptides)


prot_summary <- prot %>%
  group_by(R.Replicate, input, MS2_speed, positive) %>%
  summarize(pep_number = n_distinct(PEP.StrippedSequence)) %>%
  group_by(input, MS2_speed, positive) %>%
  summarize(mean_pep = mean(pep_number, na.rm = TRUE))

# Calculate FDR
fdr <- round((1.152 * prot_summary %>% filter(positive == FALSE) %>%
                pull(mean_pep)) /
               prot_summary %>% filter(positive == TRUE) %>%
               pull(mean_pep) * 100, 2)

# Calculate adjusted mean peptide count by species
prot_summary$y_label <- prot_summary$mean_pep

prot_summary$y_label[prot_summary$input %in% c("E. coli", "Yeast", "C. elegans", "Homo sapiens")] <-
  prot_summary$y_label[prot_summary$input %in% c("E. coli", "Yeast", "C. elegans", "Homo sapiens")] +
  ifelse(prot_summary$input == "E. coli", 790,
         ifelse(prot_summary$input == "Yeast", 520,
                ifelse(prot_summary$input == "C. elegans", 250, 0)))

# Remove Homo sapiens from y_label
prot_summary$y_label[prot_summary$input == "Homo sapiens"] <- NA


# Plot -------------------------------------------------------------------------

# Define the labels and y positions
label_data <- data.frame(
  label = c(
    rep("misidentified\nn = 10", 3), "Homo sapiens",
    rep("misidentified\nn = 9", 3), "Homo sapiens",
    rep("misidentified\nn = 14", 3), "Homo sapiens",
    rep("misidentified\nn = 14", 3), "Homo sapiens",
    rep("misidentified\nn = 2", 3), "Homo sapiens",
    rep("misidentified\nn = 24", 3), "Homo sapiens",
    rep("misidentified\nn = 13", 3), "Homo sapiens",
    rep("misidentified\nn = 17", 3), "Homo sapiens",
    rep("misidentified\nn = 8", 3), "Homo sapiens",
    rep("misidentified\nn = 5", 3), "Homo sapiens"
  ),
  y_pos = c(
    rep(2794, 3), NA,
    rep(3028, 3), NA,
    rep(3345, 3), NA,
    rep(2558, 3), NA,
    rep(1089, 3), NA,
    rep(5909, 3), NA,
    rep(5423, 3), NA,
    rep(5086, 3), NA,
    rep(3311, 3), NA,
    rep(1287, 3), NA
  ) + 100
)

# Define the FDR labels and positions
fdr_labels <- paste0(sprintf("%.2f", fdr), "%")
fdr_label_pos <- rep(c(6900, NA, NA, NA), 10)

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



