# -- Figure 4B -----------------------------------------------------------------
# The effect of library generation approach on LIT-DIA dataset characteristics.
#  The scatter in individual fragment ion peak height reproducibility between
# two replicate LIT-DDA OHPF libraries from 50 ng of starting material.

ms1 <- read.delim2("../data/20221129_073506_SN17_TestQuant_NoNormalisation_LITDIA_1_10_100_HeLa_MS1_Quant_Report.tsv")
ms2 <- read.delim2("../data/20221129_073622_SN17_TestQuant_NoNormalisation_LITDIA_1_10_100_HeLa_MS2_Quant_Report.tsv")
ms1$ms <- "MS1"
ms2$ms <- "MS2"

df<- rbind(ms1, ms2) %>%
  group_by(R.Replicate, PEP.StrippedSequence) %>%
  summarise(mean_peakh = mean(F.PeakHeight, na.rm = TRUE)) %>%
  pivot_wider(names_from = R.Replicate,
              values_from = mean_peakh)

figure_4b <- df %>%
  ggplot(aes(x = log10(`1`), y = log10(`2`)))+
  geom_pointdensity(size = 1.5)+
  scale_color_viridis_c()+
  stat_cor(method = "spearman",
           aes(label = ..r.label..),
           r.digits = 3,
           size = 5,
           label.x.npc = 0,
           label.y.npc = 0.9) +
  labs(x = "log10(F.PeakHeight) Replicate #1",
       y = "log10(F.PeakHeight) Replicate #2")+
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11, vjust = 1.2),
        strip.background = element_rect(fill = "#E5E4F7"),
        plot.title = element_text(size = 11,
                                  hjust = -0.32,
                                  vjust = 3,
                                  face = "bold"
        )
  )
