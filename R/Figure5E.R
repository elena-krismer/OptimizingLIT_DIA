# Figure 5E
# https://uclouvain-cbio.github.io/SCP.replication/articles/SCoPE2.html

# correlation between 12 single cells and bulk data
library(tidyverse)

directDIA_bulk <- read.delim2("data/figure5/20221229_095106_SN17_directDIA_LITDIA_scMS_10cells_Report.tsv")
directDIA_single<- read.delim2("data/figure5/20221229_093524_SN17_directDIA_LITDIA_scMS_12runs_Report.tsv")

directDIA_bulk["replicate"]  <-  "bulk"

directDIA_bulk <- directDIA_bulk %>%
  select(F.PeakHeight, PEP.StrippedSequence, replicate) #%>%

directDIA_single <- directDIA_single %>%
  mutate(replicate = case_when(
    grepl("23$",  R.FileName) ~ "01",
    grepl("21$",  R.FileName) ~ "02",
    grepl("19$",  R.FileName) ~ "03",
    grepl("18$",  R.FileName) ~ "04",
    grepl("17$",  R.FileName) ~ "05",
    grepl("14$",  R.FileName) ~ "06",
    grepl("13$",  R.FileName) ~ "07",
    grepl("12$",  R.FileName) ~ "08",
    grepl("07$",  R.FileName) ~ "09",
    grepl("06$",  R.FileName) ~ "10",
    grepl("05$",  R.FileName) ~ "11",
    grepl("02$",  R.FileName) ~ "12",
    TRUE ~ NA
  ))  %>%
  dplyr::select(c(F.PeakHeight, PEP.StrippedSequence, replicate))


# -- plot ----------------------------------------------------------------------
library(ggpointdensity)
library(ggpubr)

directDIA_bulk <- directDIA_bulk %>%
  dplyr::group_by(PEP.StrippedSequence) %>%
  dplyr::summarise(across(c(F.PeakHeight), median))

df5e <- directDIA_single %>%
  dplyr::group_by(PEP.StrippedSequence) %>%
  dplyr::summarise(across(c(F.PeakHeight), median)) %>%
  mutate(single = F.PeakHeight) %>%
  select(!F.PeakHeight) %>%
  merge(directDIA_bulk %>% dplyr::select(PEP.StrippedSequence, F.PeakHeight)) %>%
  mutate(bulk = F.PeakHeight)


figure_5e <- df5e%>%
  drop_na() %>%
  ggplot(aes(x = log2(bulk), y = log2(single)))+
  geom_pointdensity(size = 2.5)+
  scale_color_viridis_c()+
  stat_cor(method = "pearson",
           aes(label = ..r.label..),
           r.digits = 3,
           size = 5,
           label.x.npc = 0,
           label.y.npc = 0.9) +
  xlab(bquote("Bulk log[2]"))+
  ylab(bquote("Single cell log[2]"))+
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


# -supplement ------------------------------------------------------------------
# correlation of 12 single cells

library(ggpointdensity)
library(ggpubr)
library(GGally)

df_supplement <- directDIA_single %>%
  dplyr::group_by(PEP.StrippedSequence, replicate) %>%
  dplyr::summarise(across(c(F.PeakHeight), median)) %>%
  ungroup() %>%
  pivot_wider(names_from=replicate, values_from = F.PeakHeight)


my_fn_ggpointdensity <- function(data, mapping, N=100, ...){

  p <- ggplot(data, mapping) +
    geom_pointdensity(size = 1.5)+
    scale_color_viridis_c()
  p
}

figure_supplement <- df_supplement %>%
  dplyr::select(!PEP.StrippedSequence) %>%
  log2() %>%
  ggpairs(lower=list(continuous=my_fn_ggpointdensity )) +
  theme_bw()


