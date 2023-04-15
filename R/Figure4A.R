# -- Figure 4 A ----------------------------------------------------------------
#  The effect of library generation approach on LIT-DIA dataset characteristics
# The number of detections from 1, 10, and 100 ng of HeLa protein lysate
# injections based on library generation approach. The libraries were generated
# by either GPF LIT-DIA using 40 ng of starting material, or OHPF LIT-DDA using
# 50, 100, or 1000 ng of starting material. Error bars at the height of the bars
# indicate the spread of values in the triplicate, while error bars below the
# bars indicate the spread of high-confidence peptide detections with intensity
# coefficient of variances below 20%.

# load data --------------------------------------------------------------------

DDA_lib50_1 <- read.delim2(paste0(path, "20221213_155402_SN17_LIB_DDA_50ng_LITDIA_MS1_Enhanced_MS2_Normal_HeLa_01ng_Report.tsv"))
DDA_lib50_10 <- read.delim2(paste0(path, "20221213_154425_SN17_LIB_DDA_50ng_LITDIA_MS1_Normal_MS2_Rapid_HeLa_10ng_Report.tsv"))
DDA_lib50_100 <- read.delim2(paste0(path, "20221213_153430_SN17_LIB_DDA_50ng_LITDIA_MS1_Rapid_MS2_Rapid_HeLa_100ng_Report.tsv"))

DDA_lib100_1 <- read.delim2(paste0(path, "20221213_155257_SN17_LIB_DDA_100ng_LITDIA_MS1_Enhanced_MS2_Normal_HeLa_01ng_Report.tsv"))
DDA_lib100_10 <- read.delim2(paste0(path, "20221213_154254_SN17_LIB_DDA_100ng_LITDIA_MS1_Normal_MS2_Rapid_HeLa_10ng_Report.tsv"))
DDA_lib100_100 <- read.delim2(paste0(path, "20221213_153342_SN17_LIB_DDA_100ng_LITDIA_MS1_Rapid_MS2_Rapid_HeLa_100ng_Report.tsv"))

DDA_lib1000_1 <- read.delim2(paste0(path, "20221213_155201_SN17_LIB_DDA_1000ng_LITDIA_MS1_Enhanced_MS2_Normal_HeLa_01ng_Report.tsv"))
DDA_lib1000_10 <- read.delim2(paste0(path, "20221213_154154_SN17_LIB_DDA_1000ng_LITDIA_MS1_Normal_MS2_Rapid_HeLa_10ng_Report.tsv"))
DDA_lib1000_100 <- read.delim2(paste0(path, "20221213_153255_SN17_LIB_DDA_1000ng_LITDIA_MS1_Rapid_MS2_Rapid_HeLa_100ng_Report.tsv"))

##GFP
GFP_lib_1 <- read.delim2(paste0(path, "20221213_155119_SN17_LIB_GPF_LITDIA_MS1_Enhanced_MS2_Normal_HeLa_01ng_Report.tsv"))
GFP_lib_10 <- read.delim2(paste0(path, "20221213_154518_SN17_LIB_GPF_LITDIA_MS1_Normal_MS2_Rapid_HeLa_10ng_Report.tsv"))
GFP_lib_100 <- read.delim2(paste0(path, "20221213_153518_SN17_LIB_GPF_10ng_LITDIA_MS1_Rapid_MS2_Rapid_HeLa_100ng_Report.tsv"))


# data preparation -------------------------------------------------------------

DDA_lib50_1$fractionation <- DDA_lib50_10$fractionation <- DDA_lib50_100$fractionation <-
  DDA_lib100_1$fractionation <- DDA_lib100_10$fractionation <- DDA_lib100_100$fractionation <-
  DDA_lib1000_1$fractionation <- DDA_lib1000_10$fractionation <- DDA_lib1000_100$fractionation <- "OHPF"

GFP_lib_1$fractionation <- GFP_lib_10$fractionation <- GFP_lib_100$fractionation <- "GPF"

#library input
DDA_lib50_1$lib_input <- DDA_lib50_10$lib_input <- DDA_lib50_100$lib_input <- 50
DDA_lib100_1$lib_input <- DDA_lib100_10$lib_input <- DDA_lib100_100$lib_input <- 100
DDA_lib1000_1$lib_input <- DDA_lib1000_10$lib_input <- DDA_lib1000_100$lib_input <- 1000
GFP_lib_1$lib_input <- GFP_lib_10$lib_input <- GFP_lib_100$lib_input <- 40

#input
DDA_lib50_1$input <- DDA_lib100_1$input <- DDA_lib1000_1$input <- GFP_lib_1$input <- 1
DDA_lib50_10$input <- DDA_lib100_10$input <- DDA_lib1000_10$input <- GFP_lib_10$input <- 10
DDA_lib50_100$input <- DDA_lib100_100$input <- DDA_lib1000_100$input <- GFP_lib_100$input <- 100

lib_based <-
  rbind(DDA_lib50_1, DDA_lib50_10, DDA_lib50_100,
        DDA_lib100_1, DDA_lib100_10, DDA_lib100_100,
        DDA_lib1000_1, DDA_lib1000_10, DDA_lib1000_100,
        GFP_lib_1, GFP_lib_10, GFP_lib_100)

#GPF
GFP_lib_1_r <- GFP_lib_1 %>%
  calculate_errorbar()
GFP_lib_1_r_cv <- GFP_lib_1 %>%
  calculate_errorbar_cv()

GFP_lib_10_r <- GFP_lib_10 %>%
  calculate_errorbar()
GFP_lib_10_r_cv <- GFP_lib_10 %>%
  calculate_errorbar_cv()

GFP_lib_100_r <- GFP_lib_100 %>%
  calculate_errorbar()
GFP_lib_100_r_cv <- GFP_lib_100 %>%
  calculate_errorbar_cv()

# OHPF 50
DDA_lib50_1_r <- DDA_lib50_1 %>%
  calculate_errorbar()
DDA_lib50_1_r_cv <- DDA_lib50_1 %>%
  calculate_errorbar_cv()

DDA_lib50_10_r <- DDA_lib50_10 %>%
  calculate_errorbar()
DDA_lib50_10_r_cv <- DDA_lib50_10 %>%
  calculate_errorbar_cv()

DDA_lib50_100_r <- DDA_lib50_100 %>%
  calculate_errorbar()
DDA_lib50_100_r_cv <- DDA_lib50_100 %>%
  calculate_errorbar_cv()

# OHPF 100
DDA_lib100_1_r <- DDA_lib100_1 %>%
  calculate_errorbar()
DDA_lib100_1_r_cv <- DDA_lib100_1 %>%
  calculate_errorbar_cv()

DDA_lib100_10_r <- DDA_lib100_10 %>%
  calculate_errorbar()
DDA_lib100_10_r_cv <- DDA_lib100_10 %>%
  calculate_errorbar_cv()

DDA_lib100_100_r <- DDA_lib100_100 %>%
  calculate_errorbar()
DDA_lib100_100_r_cv <- DDA_lib100_100 %>%
  calculate_errorbar_cv()

# OHPF 1000
DDA_lib1000_1_r <- DDA_lib1000_1 %>%
  calculate_errorbar()
DDA_lib1000_1_r_cv <- DDA_lib1000_1 %>%
  calculate_errorbar_cv()

DDA_lib1000_10_r <- DDA_lib1000_10 %>%
  calculate_errorbar()
DDA_lib1000_10_r_cv <- DDA_lib1000_10 %>%
  calculate_errorbar_cv()

DDA_lib1000_100_r <- DDA_lib1000_100 %>%
  calculate_errorbar()
DDA_lib1000_100_r_cv <- DDA_lib1000_100 %>%
  calculate_errorbar_cv()


lib <- c(40, 40, 40, 50,50,50, 100,100,100, 1000, 1000, 1000)
input <- c(1,10, 100, 1, 10, 100, 1,10,100, 1,10,100)
mean <- c(
  GFP_lib_1_r["mean"],
  GFP_lib_10_r["mean"],
  GFP_lib_100_r["mean"],
  DDA_lib50_1_r["mean"],
  DDA_lib50_10_r["mean"],
  DDA_lib50_100_r["mean"],
  DDA_lib100_1_r["mean"],
  DDA_lib100_10_r["mean"],
  DDA_lib100_100_r["mean"],
  DDA_lib1000_1_r["mean"],
  DDA_lib1000_10_r["mean"],
  DDA_lib1000_100_r["mean"]
) %>% unlist()
mean_cv <- c(
  GFP_lib_1_r_cv["mean"],
  GFP_lib_10_r_cv["mean"],
  GFP_lib_100_r_cv["mean"],
  DDA_lib50_1_r_cv["mean"],
  DDA_lib50_10_r_cv["mean"],
  DDA_lib50_100_r_cv["mean"],
  DDA_lib100_1_r_cv["mean"],
  DDA_lib100_10_r_cv["mean"],
  DDA_lib100_100_r_cv["mean"],
  DDA_lib1000_1_r_cv["mean"],
  DDA_lib1000_10_r_cv["mean"],
  DDA_lib1000_100_r_cv["mean"]
) %>% unlist()
sd <- c(
  GFP_lib_1_r["sd"],
  GFP_lib_10_r["sd"],
  GFP_lib_100_r["sd"],
  DDA_lib50_1_r["sd"],
  DDA_lib50_10_r["sd"],
  DDA_lib50_100_r["sd"],
  DDA_lib100_1_r["sd"],
  DDA_lib100_10_r["sd"],
  DDA_lib100_100_r["sd"],
  DDA_lib1000_1_r["sd"],
  DDA_lib1000_10_r["sd"],
  DDA_lib1000_100_r["sd"]
) %>% unlist()
sd_cv <- c(
  GFP_lib_1_r_cv["sd"],
  GFP_lib_10_r_cv["sd"],
  GFP_lib_100_r_cv["sd"],
  DDA_lib50_1_r_cv["sd"],
  DDA_lib50_10_r_cv["sd"],
  DDA_lib50_100_r_cv["sd"],
  DDA_lib100_1_r_cv["sd"],
  DDA_lib100_10_r_cv["sd"],
  DDA_lib100_100_r_cv["sd"],
  DDA_lib1000_1_r_cv["sd"],
  DDA_lib1000_10_r_cv["sd"],
  DDA_lib1000_100_r_cv["sd"]
) %>% unlist()

df_lib_based <- data.frame(
  input = input,
  lib = lib,
  sd = sd,
  mean = mean,
  sd_cv = sd_cv,
  mean_cv = mean_cv
)

# plot -------------------------------------------------------------------------

figure_4a <- df_lib_based  %>%
  arrange(-lib) %>%
  ggplot(aes(x = as.factor(input), y = mean, fill=as.factor(lib)))+
  geom_bar(position="dodge", stat = "identity", width = 0.7, alpha=0.7, color = "black") +
  labs(y = "Identified peptides",
       x = "Input (ng)",
       fill = "Library input (ng)") +
  theme(axis.text = element_text(size = axis_text),
        axis.title = element_text(size = axis_title),
        strip.text = element_text(size = axis_text),
        strip.background = element_rect(fill = "#E5E4F7"),
        legend.text = element_text(size = axis_text),
        legend.title = element_text(size= axis_title),
        legend.position = "top")+
  scale_fill_manual(values = c("#772173", "#77a8a8", "#009599", "#005358")) +
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd),
                width = 0.2,
                size = 1, position=position_dodge(.7)) +
  geom_errorbar(aes(ymin = mean_cv - sd_cv,
                    ymax = mean_cv + sd_cv),
                width = 0.2,
                size = 1, position=position_dodge(.7))






