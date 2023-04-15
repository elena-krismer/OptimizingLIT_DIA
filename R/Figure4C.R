# ----------- Figure 4 c -------------------------------------------------------
# The effect of library generation approach on LIT-DIA dataset characteristics.
# The count of detected peptides sorted by the total number of fractions required
# to assign at least 75% of the total peptide intensity.

# load data --------------------------------------------------------------------

lib_50 <- read_tsv(paste0(path, "figure3c/evidence_50ng.txt"))
lib_100 <- read_tsv(paste0(path, "figure3c/evidence_100ng.txt"))
lib_1000 <- read_tsv(paste0(path, "figure3c/evidence_1000ng.txt"))

lib_50$input <- "50 ng"
lib_100$input <- "100 ng"
lib_1000$input <- "1000 ng"

lib <- rbind(lib_50, lib_100, lib_1000)

lib$input <- factor(lib$input, levels = c("50 ng", "100 ng", "1000 ng"))

colnames(lib) <- gsub(" ", "_", colnames(lib))

lib$Fraction <- lib$Raw_file
lib$Fraction <- sub("_LITDDA.*", "", lib$Fraction)
lib$Fraction <- sub(".*60sec_0", "", lib$Fraction)

### Count total peptides for each input
total <- lib %>%
  group_by(Sequence, input) %>%
  summarise(n_fract = length(unique(Fraction))) %>%
  group_by(input) %>%
  count(n_fract) %>%
  group_by(input) %>%
  summarise(tot = sum(n)) %>%
  pull(tot)

total <- c(rep(total[1], 8), rep(total[2], 8), rep(total[3], 8))


### Count unique peptide in each fractiion
lib_fract_count <- lib %>%
  group_by(Sequence, input) %>%
  summarise(n_fract = length(unique(Fraction))) %>%
  group_by(input) %>%
  count(n_fract) %>%
  as.data.frame()


### Calculate protortion for fractions
proportion <- lib %>%
  group_by(Sequence, input) %>%
  summarise(n_fract = length(unique(Fraction))) %>%
  group_by(input) %>%
  count(n_fract) %>%
  as.data.frame() %>%
  mutate(prop = n/total) %>%
  pull(prop)

# plot -------------------------------------------------------------------------
figure_4c <- lib_fract_count %>%
  ggplot(aes(x = n_fract, y = n, fill = input))+
  geom_col(alpha=0.7, color = "black")+
  facet_wrap(~input)+
  geom_text(aes(label = paste0(round(proportion*100, 1), "%")),
            vjust = -0.5,
            size = 3.5) +
  labs(y = "Identified peptides",
       x = "Fraction distribution")+
  theme(axis.title = element_text(size = axis_title),
        axis.text = element_text(size = axis_text),
        strip.text.x = element_text(size = axis_text),
        strip.text.y = element_text(size = axis_text, vjust = 1.2),
        strip.background = element_rect(fill = "#E5E4F7"),
        legend.position = "none")+
  scale_fill_manual(values = c("#77a8a8", "#009599", "#005358"))+
  scale_x_continuous(breaks = c(1:8))
