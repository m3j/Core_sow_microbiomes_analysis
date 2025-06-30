library(phyloseq)

source("code/Functions/stats_utils.R")
source("code/Functions/dist_utils.R")

# get phyloseq object
load("data/physeq_objects_filt0.1.RData")

# create robust aitchison distance matrix
distance_matix <- robust_ait(ps)

# define features to run permanova on
features <- c("sampling_stage", "Country", "runid", "feed_combined", "probiotic_any_binary", "antibiotic_binary", "parity_index")

# define empty dataframe to contain the results
perm_res <- data.frame()

# run permanova on each feature and add results as a dataframe
for (fill in features) {
    res <- get_permanova(physeq_obj = ps, fill = fill, dist = distance_matix)
    res_model_residual <- res[1:2, ]
    res_fill <- bind_rows(data.frame(fill), res_model_residual)
    res_fill$SumOfSqs <- round(res_fill$SumOfSqs, 0)
    res_fill$R2 <- round(res_fill$R2, 3)
    res_fill$F <- round(res_fill$F, 3)
    res_fill$"Pr(>F)" <- round(res_fill$"Pr(>F)", 4)

    # add result to data frame
    perm_res <- bind_rows(perm_res, res_fill[1:6])
}

# save results as tsv file
write_tsv(perm_res, "tables/table_S3.tsv")
