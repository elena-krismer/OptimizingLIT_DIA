# Figure 5E
# https://uclouvain-cbio.github.io/SCP.replication/articles/SCoPE2.html

# correlation between 12 single cells and bulk data


directDIA_bulk <- read.delim2("data/figure5/20221229_095106_SN17_directDIA_LITDIA_scMS_10cells_Report.tsv")
dia_40_bulk <- read.delim2("data/figure5/20221229_095124_SN17_DIA_GPF_10cells_LITDIA_scMS_10cells_Report.tsv")
dia_46_bulk <- read.delim2("data/figure5/20221229_095139_SN17_DIA_GPF_16cells_LITDIA_scMS_10cells_Report.tsv")
dia_1_bulk <- read.delim2("data/figure5/20221229_095151_SN17_DIA_GPF_10ng_LITDIA_scMS_10cells_Report.tsv")
dda_50_bulk <- read.delim2("data/figure5/20221229_095209_SN17_DDA_OHPF_50ng_LITDIA_scMS_10cells_Report.tsv")
dda_100_bulk <- read.delim2("data/figure5/20221229_095228_SN17_DDA_OHPF_100ng_LITDIA_scMS_10cells_Report.tsv")
dda_1000_bulk <- read.delim2("data/figure5/20221229_095246_SN17_DDA_OHPF_1000ng_LITDIA_scMS_10cells_Report.tsv")

directDIA_single <- read.delim2("data/figure5/20221229_093949_SN17_directDIA_LITDIA_scMS_3runs_Report.tsv")
dia_40_single <- read.delim2("data/figure5/20221229_094002_SN17_GPF_DIA_10cells_LITDIA_scMS_1cell_Report.tsv")
dia_46_single <- read.delim2("data/figure5/20221229_094017_SN17_GPF_DIA_16cells_LITDIA_scMS_1cell_Report.tsv")
dia_1_single <- read.delim2("data/figure5/20221229_094032_SN17_GPF_DIA_10ng_LITDIA_scMS_1cell_Report.tsv")
dda_50_single <- read.delim2("data/figure5/20221229_094048_SN17_OHPF_DDA_50ng_LITDIA_scMS_1cell_Report.tsv")
dda_100_single <- read.delim2("data/figure5/20221229_094613_SN17_OHPF_DDA_100ng_LITDIA_scMS_1cell_Report.tsv")
dda_1000_single <- read.delim2("data/figure5/20221229_094126_SN17_OHPF_DDA_1000ng_LITDIA_scMS_1cell_Report.tsv")

# calculate CVs of relative precursors
# "F.PeakHeight"

# -- bulk --------------------------------------------------------------------
bulk <- list(directDIA_bulk,dia_40_bulk, dia_46_bulk,
  dia_1_bulk, dda_50_bulk, dda_100_bulk, dda_1000_bulk,
  directDIA_single,dia_40_single, dia_46_single,
  dia_1_single, dda_50_single, dda_100_single, dda_1000_single)

method_label <- c("directDIA", "dia_40", "dia_46",
                "dia_1", "dda_50", "dda_100", "dda_1000",
                "directDIA", "dia_40", "dia_46",
                "dia_1", "dda_50", "dda_100", "dda_1000" )

amount <- c("bulk", "bulk", "bulk",
          "bulk", "bulk", "bulk", "bulk",
          "single", "single", "single",
          "single", "single", "single", "single")

label <- c("directDIA_bulk", "dia_40_bulk", "dia_46_bulk",
  "dia_1_bulk", "dda_50_bulk", "dda_100_bulk", "dda_1000_bulk",
  "directDIA_single", "dia_40_single", "dia_46_single",
  "dia_1_single", "dda_50_single", "dda_100_single", "dda_1000_single"
  )

for(i in 1:length(bulk)){
  bulk[[i]]["name"] <-  label[i]
  bulk[[i]]["method"] <-  method_label[i]
  bulk[[i]]["bulk_single"] <-  amount[i]
  bulk[[i]] <- bulk[[i]] %>%
    select(name, F.PeakHeight,
           PEP.StrippedSequence,
           method,
           bulk_single
           #PG.ProteinGroups
           ) %>%
    unite("pep_method", PEP.StrippedSequence:method )
}

df <- do.call(rbind,bulk)


# -- single --------------------------------------------------------------------

# median quantitative variability for every single cell


df <- df %>%
  select(!name) %>%
  group_by(pep_method, bulk_single) %>% summarise_all(funs(median(., na.rm = TRUE))) %>%
  pivot_wider(names_from=bulk_single, values_from = F.PeakHeight)


# -- plot ----------------------------------------------------------------------
library(ggpointdensity)
library(ggpubr)

figure_5e <- df %>%
  drop_na() %>%
  ggplot(aes(x = log2(bulk), y = log2(single)))+
  geom_pointdensity(size = 1.5)+
  scale_color_viridis_c()+
  stat_cor(method = "pearson",
           aes(label = ..r.label..),
           r.digits = 3,
           size = 5,
           label.x.npc = 0,
           label.y.npc = 0.9) +
  labs(x = "bulk",
       y = "single")+
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

