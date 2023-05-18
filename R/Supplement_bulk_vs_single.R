# Figure si
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


df <- rbind(directDIA_bulk, directDIA_single)

# -- plot ----------------------------------------------------------------------
# -supplement ------------------------------------------------------------------
# correlation of 12 single cells and bulk

library(ggpointdensity)
library(ggpubr)
library(GGally)

df_supplement <- df%>%
  dplyr::group_by(PEP.StrippedSequence, replicate) %>%
  dplyr::summarise(across(c(F.PeakHeight), median)) %>%
  ungroup() %>%
  pivot_wider(names_from=replicate, values_from = F.PeakHeight)


my_fn_ggpointdensity <- function(data, mapping, N=100, ...){

  p <- ggplot(data, mapping) +
    geom_pointdensity(size = 1.5)+
    scale_color_viridis_c()+
    theme(strip.text.x = element_text(size = 6),
          strip.text.y = element_text(size = 6)) +
    scale_y_continuous(breaks=c(2,5,8), limits=c(0,11)) +
    scale_x_continuous(breaks=c(2,5,8), limits=c(0,11))
  p
}

figure_supplement <- df_supplement %>%
  dplyr::select(!PEP.StrippedSequence) %>%
  log2() %>%
  ggpairs(lower=list(continuous=my_fn_ggpointdensity ),
          upper = list(continuous = wrap("cor", size = 3, colour="black")), colour="black") +
  theme(strip.text.x = element_text(size = 5),
        strip.text.y = element_text(size = 5)) +
  theme_bw()

ggsave("figure6_supplement.png", figure_supplement, height = 9, width = 9, dpi = 320)
