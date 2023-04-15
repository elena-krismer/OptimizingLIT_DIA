# - Figure 2C ------------------------------------------------------------------
# --- Repeatability of LIT-DIA measurements at different loads. ----------------
# The overall distribution of summed peak areas at either the precursor or
# protein levels.

ms1 <- read.delim2("data/20221129_073506_SN17_TestQuant_NoNormalisation_LITDIA_1_10_100_HeLa_MS1_Quant_Report.tsv")
ms2 <- read.delim2("data/20221129_073622_SN17_TestQuant_NoNormalisation_LITDIA_1_10_100_HeLa_MS2_Quant_Report.tsv")

ms1$ms <- "MS1"
ms2$ms <- "MS2"

keep_columns <- c("R.FileName", "FG.Quantity", "F.PeakArea", "ms", "PG.Quantity")

ms1_ms2 <- rbind(dplyr::select(ms1, keep_columns), dplyr::select(ms2,keep_columns))

ms1_ms2 <- ms1_ms2 %>%
  mutate(ng = case_when(
    grepl("01ng", R.FileName) ~ "1ng",
    grepl("10ng", R.FileName) ~ "10ng",
    grepl("100ng", R.FileName) ~ "100ng",
  ))


# plot -------------------------------------------------------------------------

pg_quantity <- ms1_ms2 %>%
  ggplot(aes(x=reorder(ng, PG.Quantity, na.rm = TRUE), y=PG.Quantity, fill=ms)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(trans='log10') +
  labs(y="log10(PG.Quantity)") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = axis_text),
    axis.title = element_text(size = axis_title),
    strip.text = element_text(size = axis_text),
    strip.background = element_rect(fill = "#E5E4F7"),

  ) +
  scale_fill_manual(values = c("#009599",  "#005358")) +
  facet_wrap(~ ms) +
  facet_grid()

fg_quantity <- ms1_ms2 %>%
  ggplot(aes(x=reorder(ng, FG.Quantity, na.rm = TRUE), y=PG.Quantity, fill=ms)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(trans='log10') +
  labs(y="log10(PG.Quantity)") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = axis_text),
    axis.title = element_text(size = axis_title),
    strip.text = element_text(size = axis_text),
    strip.background = element_rect(fill = "#E5E4F7"),

  ) +
  scale_fill_manual(values = c("#009599",  "#005358")) +
  facet_wrap(~ ms) +
  facet_grid()
