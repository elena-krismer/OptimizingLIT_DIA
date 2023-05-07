
compute_CV <- function(df){
  df_new <- df %>%
    mutate(sd = apply(df[, 2:(ncol(df)-2)],
                      FUN = sd, MARGIN =  1, na.rm = TRUE),
           CV = sd/mean_quantity)
  colnames(df_new)[2:(ncol(df)-2)] <- paste0("S", 1:(ncol(df)-3))
  return(df_new)
}


calculate_errorbar <- function(df){
  s1 <- df %>% select(S1) %>% drop_na() %>% nrow()
  s2 <- df %>% select(S2) %>% drop_na() %>% nrow()
  s3 <- df %>% select(S3) %>% drop_na() %>% nrow()
  all <- c(s1, s2,s3)
  result <- list()
  result["sd"] <- sd(all)
  result["mean"] <- mean(all)
  return(result)
}

calculate_errorbar_cv <- function(df){
  s1 <- df %>% select(S1, CV) %>% filter(CV < 0.2) %>% drop_na() %>% nrow()
  s2 <- df %>% select(S2, CV) %>% filter(CV < 0.2) %>% drop_na() %>% nrow()
  s3 <- df %>% select(S3, CV) %>% filter(CV < 0.2) %>% drop_na() %>% nrow()
  all <- c(s1, s2,s3)
  result <- list()
  result["sd"]  <- sd(all)
  result["mean"] <- mean(all)
  return(result)
}









filter_pep <- function(x) {
  x$FG.Quantity[as.logical(x$EG.IsImputed)] <- NA
  x_clean_with_mean <- x %>%
    group_by(R.FileName, PEP.StrippedSequence) %>%
    summarise(mean_quantity = mean(FG.Quantity, na.rm = TRUE)) %>%
    pivot_wider(names_from = R.FileName,
                values_from = mean_quantity) %>%
    suppressMessages()
  x_filter_highlight <- x_clean_with_mean %>%
    mutate(na_nbr = rowSums(is.na(x_clean_with_mean[-1]))) %>%
    mutate(mean_quantity =
             rowMeans(x_clean_with_mean[2:ncol(x_clean_with_mean)],
                      na.rm = TRUE)) %>%
    mutate(filtered = !na_nbr <= 1) %>%
    select(-na_nbr) %>%
    suppressMessages("`summarise()` has grouped output by 'R.FileName'. You can override using the `.groups` argument.")
  message("
Highlighted ", sum(x_filter_highlight$filtered) , " peptide(s) found in less than ", ncol(x_clean_with_mean) -1, " replicates")
  return(x_filter_highlight)

}

filter_protein_group <- function(x) {
  x$FG.Quantity[as.logical(x$EG.IsImputed)] <- NA
  x_clean_with_mean <- x %>%
    group_by(R.FileName, PG.ProteinGroups) %>%
    summarise(mean_quantity = mean(FG.Quantity, na.rm = TRUE)) %>%
    pivot_wider(names_from = R.FileName,
                values_from = mean_quantity) %>%
    suppressMessages()
  x_filter_highlight <- x_clean_with_mean %>%
    mutate(na_nbr = rowSums(is.na(x_clean_with_mean[-1]))) %>%
    mutate(mean_quantity =
             rowMeans(x_clean_with_mean[2:ncol(x_clean_with_mean)],
                      na.rm = TRUE)) %>%
    mutate(filtered = !na_nbr <= 1) %>%
    select(-na_nbr) %>%
    suppressMessages("`summarise()` has grouped output by 'R.FileName'. You can override using the `.groups` argument.")
  message("
Highlighted ", sum(x_filter_highlight$filtered) , " peptide(s) found in less than ", ncol(x_clean_with_mean) -1, " replicates")
  return(x_filter_highlight)

}


