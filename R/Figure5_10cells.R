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

# load data --------------------------------------------------------------------
# 10 cell pools

directDIA <- read.delim2("data/figure5/20221229_095106_SN17_directDIA_LITDIA_scMS_10cells_Report.tsv")
dia_40 <- read.delim2("data/figure5/20221229_095124_SN17_DIA_GPF_10cells_LITDIA_scMS_10cells_Report.tsv")
dia_46 <- read.delim2("data/figure5/20221229_095139_SN17_DIA_GPF_16cells_LITDIA_scMS_10cells_Report.tsv")
dia_1 <- read.delim2("data/figure5/20221229_095151_SN17_DIA_GPF_10ng_LITDIA_scMS_10cells_Report.tsv")
dda_50 <- read.delim2("data/figure5/20221229_095209_SN17_DDA_OHPF_50ng_LITDIA_scMS_10cells_Report.tsv")
dda_100 <- read.delim2("data/figure5/20221229_095228_SN17_DDA_OHPF_100ng_LITDIA_scMS_10cells_Report.tsv")
dda_1000 <- read.delim2("data/figure5/20221229_095246_SN17_DDA_OHPF_1000ng_LITDIA_scMS_10cells_Report.tsv")
