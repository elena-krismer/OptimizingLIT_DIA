# - Figure 4 D -----------------------------------------------------------------
# The effect of library generation approach on LIT-DIA dataset characteristics -
# The total amount of time required to generate DIA- and DDA-based libraries,
# including fractionation and mass spectrometry analysis.

df <- data.frame(method = c("DDA-based", "DDA-based", "DIA-based"),
                 time_c = c( "OHPF", "MS analysis", "MS analysis"),
                 time = c(100, 248, 124),
                 ms_time = c(NA, "8*31 min", "4*31 min"))

figure_4d <- df %>%
  ggplot(aes(fill = factor(time_c, levels = c("OHPF", "MS analysis")),
             y = time, x = factor(method, levels = c( "DIA-based", "DDA-based" )))) +
  geom_bar(position="stack", stat="identity",
           width = 0.6, alpha = 0.7, color = "black") +
  labs(y = "Processing time (min)",
       x = "Library generation strategies",
       fill = "a") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c("#772173", "#E5E5F7")) +
  theme(axis.title = element_text(size = axis_title),
        axis.text = element_text(size = axis_text),
        legend.text = element_text(size = axis_text),
        legend.title = element_blank(),
        legend.position = "top")
