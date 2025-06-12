# identify and record outliers and values

outlier_id <- function(x, k, id){

  ub <- median(x) + k * mad(x)
  lb <- median(x) - k * mad(x)

  return(id[which(x>ub|x<lb)])

}

ocd_ytop <- read_csv(file = '~/mrs_data/ocd_ytop_processed.csv', show_col_types = FALSE)

outliers_all <- list('> 5 mad' = unlist(lapply(select(ocd_ytop, !c(id, female, age, ocd, eos)),
                                               outlier_id, k=5, id = ocd_ytop$id)),
                     '> 4 mad' = unlist(lapply(select(ocd_ytop, !c(id, female, age, ocd, eos)),
                                               outlier_id, k=4, id = ocd_ytop$id)),
                     '> 3.5 mad' = unlist(lapply(select(ocd_ytop, !c(, female, age, ocd, eos)),
                                                 outlier_id, k=3.5, id = ocd_ytop$id)),
                     '> 3 mad' = unlist(lapply(select(ocd_ytop, !c(id, female, age, ocd, eos)),
                                                             outlier_id, k=3, id = ocd_ytop$id))) %>%

  lapply(enframe, name = "region", value = "id") %>%

  enframe(name = 'distance') %>%

  unnest('value') %>%

  mutate(region = str_remove(region, pattern = '[1234567]'),
         mri = NA,
         median = NA,
         mad = NA)

for(i in 1:nrow(outliers_all)){
  outliers_all$mri[i] <-
    select(
      filter(
        ocd_ytop, id == outliers_all$id[i]),
      outliers_all$region[i])

  outliers_all$median[i] <- sapply(select(ocd_ytop, outliers_all$region[i]), median)
  outliers_all$mad[i] <- sapply(select(ocd_ytop, outliers_all$region[i]), mad)
  outliers_all$mri <- as.double(outliers_all$mri)}

outliers_all[,4:6] <- apply(outliers_all[,4:6], 2, round, digits = 3)

write.csv(outliers_all, file = '~/mrs_data/mrs_outliers_mad.csv')
