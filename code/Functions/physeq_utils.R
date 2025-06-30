library(tidyverse)
library(phyloseq)
library(microViz)
library(hash)

# update the uncultered species with taxonomy found in EzBios db
update_species <- function(original_tax_table) {
    # get EzBio results
    ez_tax_df <- readr::read_tsv("data/EzBio_species_adjusted_otu.tsv")

    # prepare by splitting tax into levels and filter away where Genus is NA
    ez_species_df <- ez_tax_df %>%
        select("OTU", "Taxonomy") %>%
        separate(Taxonomy,
            into = c(
                "Domain", "Phylum", "Class",
                "Order", "Family", "Genus", "Species"
            ),
            sep = ";"
        ) %>%
        filter(!is.na(Genus)) %>%
        as.data.frame(.)

    # join Ezbio results with original taxa and where Genus/Species differs take the results from EzBio
    tax_table_updated <- original_tax_table %>%
        rownames_to_column("OTU") %>%
        left_join(ez_species_df, by = "OTU", suffix = c("", "_new")) %>%
        mutate(Species = coalesce(Species_new, Species)) %>%
        mutate(Genus = coalesce(Genus_new, Genus)) %>%
        select(
            OTU, Domain, Phylum, Class,
            Order, Family, Genus, Species
        )

    # make sure that the genus/species has g__/s__ in front to make equal to originial taxa
    tax_table_updated_g <- tax_table_updated %>%
        mutate(Genus = if_else(str_detect(Genus, "g__"),
            Genus,
            paste0("g__", Genus)
        )) %>%
        mutate(Species = if_else(str_detect(Species, "s__"),
            Species,
            paste0("s__", Species)
        ))

    # make OTU as rownames and remove the OTU column
    rownames(tax_table_updated_g) <- tax_table_updated_g$OTU
    tax_table_updated <- tax_table_updated_g %>% select(-OTU)

    return(tax_table_updated)
}

# calculate abundance of input physeq object and return df with otu and corresponding abundance(reads and %) across the entire data
get_abundance <- function(ps_obj) {
    otu_count <- phyloseq::otu_table(ps_obj)
    total_reads <- sum(otu_count)
    abundance <- otu_count %>%
        rowSums(.) %>%
        sort(., decreasing = TRUE) %>%
        as.data.frame() %>%
        rownames_to_column("otu_id") %>%
        rename(abundance_reads = ".") %>%
        mutate("abundance_%" = round(abundance_reads / total_reads * 100, 2))

    return(abundance)
}

# calculate the prevalence of each otu in the input physeq object and return a df with OTU as row and one column with prevalence in percentage
get_prevalence <- function(ps_obj) {
    otu_counts <- otu_table(ps_obj)
    prevalence <- apply(otu_counts, 1, function(x) round(mean(x > 0) * 100, 2)) %>%
        sort(decreasing = TRUE) %>%
        as.data.frame() %>%
        rename(prevalence = ".")

    return(prevalence)
}

# split a phyloseq object into two based on input variable ( default sampling stage) and return as a list
split_physeq_in_two <- function(physeq_obj, feature, value1, value2) {
    # filter phyloseq object on feature given the different values
    physeq_1 <- microViz::ps_filter(physeq_obj, .data[[feature]] == value1)
    physeq_2 <- microViz::ps_filter(physeq_obj, .data[[feature]] == value2)

    # returns a list of the two phyloseq objects
    return(list(physeq_1 = physeq_1, physeq_2 = physeq_2))
}

# takes a phyloseq object and returns a dict with otu ids as key and taxonomy at a defined rank
otu_tax_dict <- function(physeq_obj, rank) {
    # get taxa and turn rownames into a column named otu_id
    tax <- data.frame(tax_table(physeq_obj)) %>%
        tibble::rownames_to_column(., "otu_id")
    # make column names to lowcase
    colnames(tax) <- tolower(colnames(tax))

    # instanciate an emtpy hash
    otu_hash <- hash()

    # occupy hash with otu id as key and taxa as value
    for (row in seq_len(nrow(tax))) {
        key <- tax[row, ]$otu_id
        value <- tax[row, ][[tolower(rank)]]
        otu_hash[key] <- value
    }

    return(otu_hash)
}

# will take a dataframe column with OTU ids in format OUT239 and
# return the dataframe with additional column named tax_names containing taxa at defined rank.
otu_id2tax_names_col <- function(otu_df, otu_col, physeq_obj, rank) {
    # create otu hash(key = otu, value = taxa at rank)
    otu_hash <- otu_tax_dict(physeq_obj, rank)

    # set tax names to be otu ids
    tax_names <- otu_df[otu_col]

    # loop through the length of input df
    for (idx in seq_len(length(otu_df[[otu_col]]))) {
        # get taxa from otu hash looking up using otuid as key
        tax_name <- otu_hash[[otu_df[otu_col][idx, ]]]
        # change tax_names to the taxname at index
        tax_names[idx, ] <- tax_name
    }

    # create a new column named tax_names and setting it equal to tax_names.
    otu_df["tax_names"] <- tax_names

    # return full df with additional tax column
    return(otu_df)
}

# takes one otu id and returns tax according to specified rank
otu2tax <- function(physeq_obj, rank, otu) {
    # create otu hash(key = otu, value = taxa at rank)
    otu_dict <- otu_tax_dict(physeq_obj, rank)

    # try to get taxa value given otu key return NA if key is unvalid.
    otu <- tryCatch(
        {
            get(otu, otu_dict)
        },
        error = function(e) {
            message("the otu number '", otu, "' not found. Returning NA.")
            NA
        }
    )
    return(otu)
}

# filt a phyloseq object based on a median threshold.
filt_on_median <- function(physeq_obj, median_threshold) {
    # get relative abundance across entire data and make rownames a column called taxon
    physeq_ra <- apply((otu_table(physeq_obj)), 2, function(x) x / sum(x)) %>%
        as.data.frame() %>%
        rownames_to_column("taxon")

    # get abundance columns removing taxon column
    abund_cols <- names(physeq_ra)[-1]

    # create additional column containing median abundance of each row/otu
    physeq_ra <- physeq_ra %>%
        mutate(median_ra = apply(select(., all_of(abund_cols)), 1, median))

    # filter on median column keep those above input threshold
    physeq_ra_filtered <- physeq_ra %>% filter(median_ra > median_threshold)

    return(physeq_ra_filtered)
}
