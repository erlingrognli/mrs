# identify and record outliers and values

# outlier_id <- function(x, k, id){
#
#   ub <- median(x) + k * mad(x)
#   lb <- median(x) - k * mad(x)
#
#   return(id[which(x>ub|x<lb)])
#
# }

# outliers_all <- list('> 5 mad' = unlist(lapply(select(ocd_ytop, !c(SubjID, female, age, ocd, eop)),
#                                                outlier_id, k=5, id = ocd_ytop$SubjID)),
#                      '> 4 mad' = unlist(lapply(select(ocd_ytop, !c(SubjID, female, age, ocd, eop)),
#                                                outlier_id, k=4, id = ocd_ytop$SubjID)),
#                      '> 3.5 mad' = unlist(lapply(select(ocd_ytop, !c(SubjID, female, age, ocd, eop)),
#                                                  outlier_id, k=3.5, id = ocd_ytop$SubjID))) %>%
#
#   lapply(enframe, name = "region", value = "id") %>%
#
#   enframe(name = 'distance') %>%
#
#   unnest('value') %>%
#
#   mutate(region = str_remove(region, pattern = '[1234567]'),
#          mri = NA,
#          median = NA)
#
# for(i in 1:nrow(outliers_all)){
#   outliers_all$mri[i] <-
#     select(
#       filter(
#         ocd_ytop, SubjID == outliers_all$id[i]),
#       outliers_all$region[i])
#
#   outliers_all$median[i] <- sapply(select(ocd_ytop, outliers_all$region[i]), median)
#   outliers_all$mri <- as.double(outliers_all$mri)}
#
# write.csv(outliers_all, file = '~/mrs_data/mrs_outliers_mad.csv')
