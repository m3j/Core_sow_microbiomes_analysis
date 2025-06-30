library(vegan)
library(phyloseq)

# calculate robust aichison distances on count table from physeq obj.
robust_ait <- function(physeq_obj) {
    df <- as.data.frame(otu_table(physeq_obj))
    dist <- vegan::vegdist(t(df), method = "robust.aitchison", upper = TRUE)
    return(dist)
}
