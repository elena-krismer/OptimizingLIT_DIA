# --------------- Figure 2A  ---------------------------------------------------
# ---- Repeatability of LIT-DIA measurements at different loads. ---------------
# The number of detected peptides from triplicate 1, 10, and 100 ng of HeLa
# protein lysate. Error bars at the height of the bars indicate the spread of
# values in the triplicate, while error bars below the bars indicate the spread
# of high-confidence peptide detections with intensity coefficient of variances
# below 20%.

# Define file paths
file_paths <- c("data/SN17_directDIA_LITDIA_MS1_Enhanced_MS2_Normal_HeLa_01ng_Report.tsv",
                "data/SN17_directDIA_LITDIA_MS1_Normal_MS2_Rapid_HeLa_10ng_MS1_Quant_Report.tsv",
                "data/SN17_directDIA_LITDIA_MS1_Rapid_MS2_Rapid_HeLa_100ng_MS1_Quant_Report.tsv")

# Define labels and input values for each file
file_labels <- c("1 ng \n MS1 Enhanced - MS Normal",
                 "10 ng \n MS1 Normal - MS Rapid",
                 "100 ng \n MS1 Rapid - MS Rapid")
input_values <- c("1 ng", "10 ng", "100 ng")

# Read and process each file
data_list <- lapply(file_paths, function(path) {
  data <- read.delim2(path)
  data <- data %>%
    filter_pep() %>%
    compute_CV()
  input <- input_values[file_paths == path]
  label_x <- file_labels[file_paths == path]
  data %>% mutate(input = input, label_x = label_x)
})

# Combine data and compute summary statistics
df <- do.call(rbind, data_list)
df <- df %>%
  arrange(desc(CV)) %>%
  mutate(
    cv_label = case_when(
      CV <= 0.2  ~ "< 20 %",
      CV > 0.2 ~ "> 20 %",
      TRUE ~ as.character(NA)
    )
  ) %>%
  dplyr::filter(!is.na(cv_label)) %>%
  group_by(input, cv_label, CV, label_x) %>%
  summarise(id_S1 = sum(S1),
            id_S2 = sum(S2),
            id_S3 = sum(S3)) %>%
  rowwise() %>%
  mutate(
    mean_n_peptides = mean(c(id_S1, id_S2, id_S3)),
    sd_id = sd(c(id_S1, id_S2, id_S3))
  ) %>%
  dplyr::group_by(input) %>%
  mutate(pos = cumsum(mean_n_peptides))

# Get data for plotting figure 2
df_plot <- get_plotting_data_fig2(df)


# Plot -------------------------------------------------------------------------

figure_2a <- df_plot %>%
  ggplot(aes(x=as.factor(label), y = mean)) +
  geom_bar(stat = "identity", width = 0.7,  alpha=0.7, fill = "#005358", color = "black") +
  geom_text(aes(label = label_2),  angle=90) +
  labs(x = "Input (ng)",
       y = "Identified peptides",
       fill = "CV") +
  theme_bw()  +
  theme(axis.text = element_text(size = axis_text),
        axis.title = element_text(size = axis_title),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 11),
        strip.background = element_rect(fill = "#E5E4F7"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.8, "cm"),
        legend.position = "top") +
  scale_fill_manual(values = c("#e6feff", "#005358"))
