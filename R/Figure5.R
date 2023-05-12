# - Figure 5 -------------------------------------------------------------------
# The number of peptides detected in triplicate (A) single-cell injections and
# (B) 10-cell pools searched using different library generation strategies.
# Cells were sorted via cellenONE and analyzed by LIT-DIA on Whisper 100 40 SPD
# and searched in directDIA mode, or searched against DDA-based spectral
# libraries and DIA-based chromatogram libraries generated from different input
# levels. Error bars at the height of the bars indicate the spread of values in
# the triplicate, while error bars below the bars indicate the spread of
# high-confidence peptide detections with intensity coefficient of variances
# below 20%.

# load data --------------------------------------------------------------------
# 10 cell pools

directDIA <- read.delim2("data/figure4b/20221229_095106_SN17_directDIA_LITDIA_scMS_10cells_Report.tsv")
dia_40 <- read.delim2("data/figure4b/20221229_095124_SN17_DIA_GPF_10cells_LITDIA_scMS_10cells_Report.tsv")
dia_46 <- read.delim2("data/figure4b/20221229_095139_SN17_DIA_GPF_16cells_LITDIA_scMS_10cells_Report.tsv")
dia_1 <- read.delim2("data/figure4b/20221229_095151_SN17_DIA_GPF_10ng_LITDIA_scMS_10cells_Report.tsv")
dda_50 <- read.delim2("data/figure4b/20221229_095209_SN17_DDA_OHPF_50ng_LITDIA_scMS_10cells_Report.tsv")
dda_100 <- read.delim2("data/figure4b/20221229_095228_SN17_DDA_OHPF_100ng_LITDIA_scMS_10cells_Report.tsv")
dda_1000 <- read.delim2("data/figure4b/20221229_095246_SN17_DDA_OHPF_1000ng_LITDIA_scMS_10cells_Report.tsv")

# single-cell injections
directDIA <- read.delim2("data/figure4a/20221229_093949_SN17_directDIA_LITDIA_scMS_3runs_Report.tsv")
dia_40 <- read.delim2("data/figure4a/20221229_094002_SN17_GPF_DIA_10cells_LITDIA_scMS_1cell_Report.tsv")
dia_46 <- read.delim2("data/figure4a/20221229_094017_SN17_GPF_DIA_16cells_LITDIA_scMS_1cell_Report.tsv")
dia_1 <- read.delim2("data/figure4a/20221229_094032_SN17_GPF_DIA_10ng_LITDIA_scMS_1cell_Report.tsv")
dda_50 <- read.delim2("data/figure4a/20221229_094048_SN17_OHPF_DDA_50ng_LITDIA_scMS_1cell_Report.tsv")
dda_100 <- read.delim2("data/figure4a/20221229_094613_SN17_OHPF_DDA_100ng_LITDIA_scMS_1cell_Report.tsv")
dda_1000 <- read.delim2("data/figure4a/20221229_094126_SN17_OHPF_DDA_1000ng_LITDIA_scMS_1cell_Report.tsv")


# prepare data -----------------------------------------------------------------
library(dplyr)
library(purrr)

# prepare data
datasets <- list(directDIA, dia_40, dia_46, dia_1, dda_50, dda_100, dda_1000)

processed_datasets <- datasets %>%
  map_dfr(function(dataset) {
    dataset %>%
      filter_pep() %>%
      compute_CV() %>%
      mutate(input = str_remove(names(dataset), "_s$"),
             label_x = "1 ng \n MS1 Enhanced - MS Normal") %>%
      calculate_errorbar() %>%
      calculate_errorbar_cv()
  })

# create df
df_plot <- tibble(
  label = c("directDIA", "DIA-based 40 cells", "DIA-based 64 cells", "DIA-based 40 ng",
            "DDA-based 50 ng", "DDA-based 100 ng", "DDA-based 1000 ng"),
  sd = processed_datasets %>% pull(sd),
  mean = processed_datasets %>% pull(mean),
  sd_cv = processed_datasets %>% pull(sd_cv),
  mean_cv = processed_datasets %>% pull(mean_cv)
)

# plot -------------------------------------------------------------------------

figure_4a <- df_plot %>%
  ggplot(aes(x=factor(label, level = label), y = mean)) +
  geom_bar(stat = "identity", width = 0.7,  alpha=0.7, fill = "#005358", color = "black") +
  labs(x = "Input (ng)",
       y = "Identified peptides",
       fill = "CV") +
  theme_bw()  +
  theme(axis.text = element_text(size = axis_text),
        axis.title = element_text(size = axis_title),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 11),
        strip.background = element_rect(fill = "#E5E4F7"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.8, "cm"),
        legend.position = "top")+
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd),
                width = 0.2,
                size = 1) +
  geom_errorbar(aes(ymin = mean_cv - sd_cv,
                    ymax = mean_cv + sd_cv),
                width = 0.2,
                size = 1)

# load data --------------------------------------------------------------------

directDIA <- read.delim2("data/figure4b/20221229_095106_SN17_directDIA_LITDIA_scMS_10cells_Report.tsv")
dia_40 <- read.delim2("data/figure4b/20221229_095124_SN17_DIA_GPF_10cells_LITDIA_scMS_10cells_Report.tsv")
dia_46 <- read.delim2("data/figure4b/20221229_095139_SN17_DIA_GPF_16cells_LITDIA_scMS_10cells_Report.tsv")
dia_1 <- read.delim2("data/figure4b/20221229_095151_SN17_DIA_GPF_10ng_LITDIA_scMS_10cells_Report.tsv")
dda_50 <- read.delim2("data/figure4b/20221229_095209_SN17_DDA_OHPF_50ng_LITDIA_scMS_10cells_Report.tsv")
dda_100 <- read.delim2("data/figure4b/20221229_095228_SN17_DDA_OHPF_100ng_LITDIA_scMS_10cells_Report.tsv")
dda_1000 <- read.delim2("data/figure4b/20221229_095246_SN17_DDA_OHPF_1000ng_LITDIA_scMS_10cells_Report.tsv")

# prepare data -----------------------------------------------------------------
library(dplyr)
library(purrr)

# prepare data
datasets <- list(directDIA, dia_40, dia_46, dia_1, dda_50, dda_100, dda_1000)

processed_datasets <- datasets %>%
  map_dfr(function(dataset) {
    dataset %>%
      filter_pep() %>%
      compute_CV() %>%
      mutate(input = str_remove(names(dataset), "_s$"),
             label_x = "1 ng \n MS1 Enhanced - MS Normal") %>%
      calculate_errorbar() %>%
      calculate_errorbar_cv()
  })

# create df
df_plot <- tibble(
  label = c("directDIA", "DIA-based 40 cells", "DIA-based 64 cells", "DIA-based 40 ng",
            "DDA-based 50 ng", "DDA-based 100 ng", "DDA-based 1000 ng"),
  sd = processed_datasets %>% pull(sd),
  mean = processed_datasets %>% pull(mean),
  sd_cv = processed_datasets %>% pull(sd_cv),
  mean_cv = processed_datasets %>% pull(mean_cv)
)

# plot -------------------------------------------------------------------------

figure_4a <- df_plot %>%
  ggplot(aes(x=factor(label, level = label), y = mean)) +
  geom_bar(stat = "identity", width = 0.7,  alpha=0.7, fill = "#005358", color = "black") +
  labs(x = "Input (ng)",
       y = "Identified peptides",
       fill = "CV") +
  theme_bw()  +
  theme(axis.text = element_text(size = axis_text),
        axis.title = element_text(size = axis_title),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 11),
        strip.background = element_rect(fill = "#E5E4F7"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.8, "cm"),
        legend.position = "top")+
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd),
                width = 0.2,
                size = 1) +
  geom_errorbar(aes(ymin = mean_cv - sd_cv,
                    ymax = mean_cv + sd_cv),
                width = 0.2,
                size = 1)
