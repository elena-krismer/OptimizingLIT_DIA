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

# Figure A and C

# load data --------------------------------------------------------------------

# single-cell injections
directDIA <- read.delim2("data/figure5/20221229_093949_SN17_directDIA_LITDIA_scMS_3runs_Report.tsv")
dia_40 <- read.delim2("data/figure5/20221229_094002_SN17_GPF_DIA_10cells_LITDIA_scMS_1cell_Report.tsv")
dia_46 <- read.delim2("data/figure5/20221229_094017_SN17_GPF_DIA_16cells_LITDIA_scMS_1cell_Report.tsv")
dia_1 <- read.delim2("data/figure5/20221229_094032_SN17_GPF_DIA_10ng_LITDIA_scMS_1cell_Report.tsv")
dda_50 <- read.delim2("data/figure5/20221229_094048_SN17_OHPF_DDA_50ng_LITDIA_scMS_1cell_Report.tsv")
dda_100 <- read.delim2("data/figure5/20221229_094613_SN17_OHPF_DDA_100ng_LITDIA_scMS_1cell_Report.tsv")
dda_1000 <- read.delim2("data/figure5/20221229_094126_SN17_OHPF_DDA_1000ng_LITDIA_scMS_1cell_Report.tsv")



filter_protein_group <- function(x) {
  x$PG.Quantity[as.logical(x$EG.IsImputed)] <- NA
  x_clean_with_mean <- x %>%
    group_by(R.FileName, PG.ProteinGroups) %>%
    summarise(mean_quantity = mean(PG.Quantity, na.rm = TRUE)) %>%
    pivot_wider(names_from = R.FileName,
                values_from = mean_quantity) %>%
    suppressMessages()
  x_filter_highlight <- x_clean_with_mean %>%
    mutate(na_nbr = rowSums(is.na(x_clean_with_mean[-1]))) %>%
    mutate(mean_quantity =
             rowMeans(x_clean_with_mean[2:ncol(x_clean_with_mean)],
                      na.rm = TRUE)) %>%
    mutate(filtered = !na_nbr <= 1) %>%
    select(-na_nbr) %>%
    suppressMessages("`summarise()` has grouped output by 'R.FileName'. You can override using the `.groups` argument.")
  message("
Highlighted ", sum(x_filter_highlight$filtered) , " peptide(s) found in less than ", ncol(x_clean_with_mean) -1, " replicates")
  return(x_filter_highlight)

}

# prepare data -----------------------------------------------------------------
directDIA_s <- directDIA%>%
  #filter_pep() %>%
  filter_protein_group() %>%
  compute_CV() %>%
  mutate(input = "directDIA",
         label_x = "1 ng \n MS1 Enhanced - MS Normal")
directDIA_errorbar <- directDIA_s %>% calculate_errorbar()
directDIA_errorbar_cv <- directDIA_s %>% calculate_errorbar_cv()

dia_40_s <- dia_40 %>%
  #filter_pep() %>%
  filter_protein_group() %>%
  compute_CV() %>%
  mutate(input = "dia40",
         label_x = "1 ng \n MS1 Enhanced - MS Normal")
dia_40_errorbar <- dia_40_s %>% calculate_errorbar()
dia_40_errorbar_cv <- dia_40_s %>% calculate_errorbar_cv()

dia_46_s <- dia_46 %>%
 # filter_pep() %>%
  filter_protein_group() %>%
  compute_CV() %>%
  mutate(input = "dia16",
         label_x = "1 ng \n MS1 Enhanced - MS Normal")
dia_46_errorbar <- dia_46_s %>% calculate_errorbar()
dia_46_errorbar_cv <- dia_46_s %>% calculate_errorbar_cv()

dia_1_s <- dia_1 %>%
 # filter_pep() %>%
  filter_protein_group() %>%
  compute_CV() %>%
  mutate(input = "dia",
         label_x = "1 ng \n MS1 Enhanced - MS Normal")
dia_1_errorbar <- dia_1_s %>% calculate_errorbar()
dia_1_errorbar_cv <- dia_1_s %>% calculate_errorbar_cv()

dda_50_s <- dda_50 %>%
  #filter_pep() %>%
  filter_protein_group() %>%
  compute_CV() %>%
  mutate(input = "dda50",
         label_x = "1 ng \n MS1 Enhanced - MS Normal")
dda_50_errorbar <- dda_50_s %>% calculate_errorbar()
dda_50_errorbar_cv <- dda_50_s %>% calculate_errorbar_cv()

dda_100_s <- dda_100 %>%
  #filter_pep() %>%
  filter_protein_group() %>%
  compute_CV() %>%
  mutate(input = "dda100",
         label_x = "1 ng \n MS1 Enhanced - MS Normal")
dda_100_errorbar <- dda_100_s %>% calculate_errorbar()
dda_100_errorbar_cv <- dda_100_s %>% calculate_errorbar_cv()

dda_1000_s <- dda_1000 %>%
  #filter_pep() %>%
  filter_protein_group() %>%
  compute_CV() %>%
  mutate(input = "dda1000",
         label_x = "1 ng \n MS1 Enhanced - MS Normal")
dda_1000_errorbar <- dda_1000_s %>% calculate_errorbar()
dda_1000_errorbar_cv <- dda_1000_s %>% calculate_errorbar_cv()

# create df --------------------------------------------------------------------
df_plot = data.frame()

label <- c("directDIA",
           "DIA-based 40 cells", "DIA-based 64 cells", "DIA-based 40 ng",
           "DDA-based 50 ng", "DDA-based 100 ng", "DDA-based 1000 ng")


sd <- c(directDIA_errorbar$sd,
        dia_40_errorbar$sd, dia_46_errorbar$sd, dia_1_errorbar$sd,
        dda_50_errorbar$sd, dda_100_errorbar$sd, dda_1000_errorbar$sd)

mean <- c(directDIA_errorbar$mean,
          dia_40_errorbar$mean, dia_46_errorbarr$mean,   dia_1_errorbar$mean,
          dda_50_errorbar$mean, dda_100_errorbar$mean, dda_1000_errorbar$mean)

sd_cv <- c(directDIA_errorbar_cv$sd,
           dia_40_errorbar_cv$sd, dia_46_errorbar_cv$sd,   dia_1_errorbar_cv$sd,
           dda_50_errorbar_cv$sd, dda_100_errorbar_cv$sd, dda_1000_errorbar_cv$sd)

mean_cv <- c(directDIA_errorbar_cv$mean,
             dia_40_errorbar_cv$mean, dia_46_errorbar_cv$mean,  dia_1_errorbar_cv$mean,
             dda_50_errorbar_cv$mean, dda_100_errorbar_cv$mean, dda_1000_errorbar_cv$mean)

df_plot <- data.frame(
  label = label,
  sd = sd,
  mean = mean,
  sd_cv = sd_cv,
  mean_cv = mean_cv
)


# plot -------------------------------------------------------------------------

figure_5a <- df_plot %>%
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

figure_5c <- df_plot %>%
  ggplot(aes(x=factor(label, level = label), y = mean)) +
  geom_bar(stat = "identity", width = 0.7,  alpha=0.7, fill = "#005358", color = "black") +
  labs(x = "Input (ng)",
       y = "Identified ProteinGroups",
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


ggarrange(figure_5c, figure_5d, labels=c("", ""), nrow=1, ncol=2)

