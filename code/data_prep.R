library(tidyverse)
library(qiime2R)
library(phyloseq)
source("code/Functions/physeq_utils.R")

# perpare metadata
metadata_intial <- readr::read_tsv("data/metadata.tsv")

# replace Invalid Number
metadata <- metadata_intial %>%
    mutate(across(
        where(~ is.character(.x) || is.factor(.x)),
        ~ gsub("Invalid Number", NA, as.character(.x))
    ))

# rename columns - removes () and yes/no and spaces
metadata <- metadata %>%
    rename_with(~ sub("\\(", "", .)) %>%
    rename_with(~ sub("\\)", "", .)) %>%
    rename_with(~ sub("_yes/no", "", .)) %>%
    rename_with(~ sub(" ", "_", .))

# convert to binary
metadata <- metadata %>%
    mutate(
        GD109_probiotic_binary =
            case_when(
                tolower(GD109_probiotic) == "no" ~ "0",
                tolower(GD109_probiotic) == "yes" ~ "1"
            )
    ) %>%
    mutate(
        LD21_probiotic_binary =
            case_when(
                tolower(LD21_probiotic) == "no" ~ "0",
                tolower(LD21_probiotic) == "yes" ~ "1"
            )
    ) %>%
    mutate(
        antibiotic_binary =
            case_when(
                antibiotic == "no" ~ "0",
                antibiotic == "yes" & farm_number == "ES_3" & sampling_stage == "G109" ~ "0", # only administered in lactation
                antibiotic == "yes" & farm_number == "ES_2" & sampling_stage == "G109" ~ "0", # only administered in lactation
                antibiotic == "yes" ~ "1"
            )
    ) %>%
    mutate(
        probiotic_any_binary =
            case_when(
                GD109_probiotic_binary == 0 & LD21_probiotic_binary == 0 ~ "0",
                LD21_probiotic_binary == 1 ~ "1",
                GD109_probiotic_binary == 1 ~ "1"
            )
    )

# categorizing data:
metadata <- metadata %>%
    mutate(
        parity_index =
            case_when(
                parity == 1 ~ "1",
                parity == 2 ~ "2",
                parity == 3 ~ "3",
                parity == 4 ~ "4",
                parity >= 5 ~ ">5"
            )
    ) %>%
    mutate(
        GD109_temp_gestation_unit_index =
            case_when(
                GD109_temp_gestation_unit < 19 ~ "<19",
                GD109_temp_gestation_unit >= 19 &
                    GD109_temp_gestation_unit <= 22 ~ "19 <= x <= 22",
                GD109_temp_gestation_unit > 22 ~ ">22"
            )
    ) %>%
    mutate(
        LD21_temp_lactation_unit_index =
            case_when(
                LD21_temp_lactation_unit <= 22 ~ "<=22",
                LD21_temp_lactation_unit > 22 &
                    LD21_temp_lactation_unit <= 25 ~ "22 < x <= 25",
                LD21_temp_lactation_unit > 25 ~ ">25"
            )
    ) %>%
    mutate(
        no_piglets_alive_index =
            case_when(
                no_piglets_alive < 13 ~ "<13",
                no_piglets_alive >= 13 &
                    no_piglets_alive <= 19 ~ "13 <= x <= 19",
                no_piglets_alive > 19 ~ ">19"
            )
    ) %>%
    mutate(
        no_piglets_weaned_per_year_sow_not_for_gilts_index =
            case_when(
                no_piglets_weaned_per_year_sow_not_for_gilts < 28.4 ~ "<28.4",
                no_piglets_weaned_per_year_sow_not_for_gilts >= 28.4 &
                    no_piglets_weaned_per_year_sow_not_for_gilts <= 37.1 ~ "28.4 <= x <= 37.1",
                no_piglets_weaned_per_year_sow_not_for_gilts > 37.1 ~ ">37.1"
            )
    ) %>%
    mutate(
        no_piglets_weaned_per_year_sow_not_for_gilts_index2 =
            case_when(
                no_piglets_weaned_per_year_sow_not_for_gilts_index == ">37.1" ~ "Good>37.1",
                no_piglets_weaned_per_year_sow_not_for_gilts_index == "28.4 <= x <= 37.1" ~ "Bad<37.1",
                no_piglets_weaned_per_year_sow_not_for_gilts_index == "<28.4" ~ "Bad<37.1"
            )
    ) %>%
    mutate(
        no_stillborn_piglets_index =
            case_when(
                as.numeric(no_stillborn_piglets) < 1 ~ "<1",
                as.numeric(no_stillborn_piglets) >= 1 &
                    as.numeric(no_stillborn_piglets) <= 2 ~ "1 <= x <= 2",
                as.numeric(no_stillborn_piglets) > 2 ~ ">2"
            )
    ) %>%
    mutate(
        enterotypes4 = case_when(
            sampling_stage == "G109" & no_piglets_weaned_per_year_sow_not_for_gilts_index2 == "Good>37.1" ~ "Gestation_good",
            sampling_stage == "G109" & no_piglets_weaned_per_year_sow_not_for_gilts_index2 == "Bad<37.1" ~ "Gestation_bad",
            sampling_stage == "L21" & no_piglets_weaned_per_year_sow_not_for_gilts_index2 == "Good>37.1" ~ "Lactation_good",
            sampling_stage == "L21" & no_piglets_weaned_per_year_sow_not_for_gilts_index2 == "Bad<37.1" ~ "Lactation_bad",
        )
    ) %>%
    mutate(diet_change = case_when(
        as.character(GD019_feed_formulation_major_ingredients) == as.character(LD21_feed_formulation_major_ingredients) ~ "no_change",
        as.character(GD019_feed_formulation_major_ingredients) != as.character(LD21_feed_formulation_major_ingredients) ~ "change",
        is.na(GD019_feed_formulation_major_ingredients) | is.na(LD21_feed_formulation_major_ingredients) ~ NA_character_
    ))

# make sample_id as row names for phyloseq generation
metadata <- metadata %>%
    column_to_rownames(var = "sampleid")

# factorize binary values for discrete coloring
metadata$probiotic_any_binary <- as.factor(metadata$probiotic_any_binary)
metadata$antibiotic_binary <- as.factor(metadata$antibiotic_binary)
metadata$diet_change <- as.factor(metadata$diet_change)

# reorder factor levels
metadata$parity_index <- factor(metadata$parity_index, levels = c("1", "2", "3", "4", ">5"))
metadata$no_piglets_weaned_per_year_sow_not_for_gilts_index <- factor(
    metadata$no_piglets_weaned_per_year_sow_not_for_gilts_index,
    levels = c(">37.1", "28.4 <= x <= 37.1", "<28.4")
)
metadata$LD21_temp_lactation_unit_index <- factor(
    metadata$LD21_temp_lactation_unit_index,
    levels = c("<=22", "22 < x <= 25", ">25")
)
metadata$GD109_temp_gestation_unit_index <- factor(
    metadata$GD109_temp_gestation_unit_index,
    levels = c("<19", "19 <= x <= 22", ">22")
)


# prepare taxonomy data
# get silva data
silva_tax <- qiime2R::read_qza("data/silva-138-ssu-nr99-tax.qza")
tax_df <- silva_tax$data

# get otu table
otu_table <- readr::read_tsv("data/otu_table.tsv")

# join otu data with silva data on silva accession number
otu_tax <- otu_table %>% left_join(tax_df, by = "Feature.ID")

# rename Clostridium sensu strico to Clostridium
otu_tax_edit <- otu_tax %>%
    mutate(Taxon = if_else(stringr::str_detect(Taxon, "g__Clostridium_sensu_stricto_1"),
        stringr::str_replace(Taxon, "g__Clostridium_sensu_stricto_1", "g__Clostridium"),
        Taxon
    ))

# edit column names removing X and replace . wiht -
otu_tax_edit_format <- otu_tax_edit %>%
    rename_with(~ stringr::str_replace(., "X", "")) %>%
    rename_with(~ stringr::str_replace_all(., "\\.", "-"))

# taking the silva taxonomy (last column) and make it its own dataframe naming the rows by OTU and row number
tax <- as.data.frame(otu_tax_edit_format[, length(otu_tax_edit_format)])
rownames(tax) <- paste0("OTU", seq_len(nrow(tax)))

# split the taxa from one string into the ranks
tax_split <- tax %>%
    separate(Taxon,
        into = c(
            "Domain", "Phylum", "Class",
            "Order", "Family", "Genus", "Species"
        ),
        sep = ";"
    )

# add species names found of the relevant uncultured species from EzBio
tax_split <- update_species(tax_split)

# rename from Clostridiales_bacterium to uncultured
# renaming mislabelled Trichuris trichiura to Bacillus
tax_split["OTU3765", "Species"] <- "s__uncultured_bacterium"
tax_split["OTU476", "Species"] <- "s__uncultured_bacterium"
tax_split["OTU763", "Class"] <- "c__Bacilli"
tax_split["OTU763", "Order"] <- "o__Bacillales"
tax_split["OTU763", "Family"] <- "f__Bacillaceae"
tax_split["OTU763", "Genus"] <- "g__Bacillus"
tax_split["OTU763", "Species"] <- "s__uncultured_bacterium"

# rename uncultered/unknown/empty levels(eg. f__) to be unique
for (t in c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    # define empty annotation and unknown
    empty_annotation <- paste(tolower(substr(t, 1, 1)), "__", sep = "")
    unknown <- paste("unknown_", tolower(t), sep = "")

    # get all values form the indexed column (t)
    column_values <- tax_split[, t]

    # find NA and empty annotations
    is_empty_or_na <- is.na(column_values) | column_values == empty_annotation

    # get the original label from NA and empty annotations
    base_labels <- column_values[is_empty_or_na]

    # add the unknown label to NA
    base_labels[is.na(base_labels)] <- unknown

    # make base_labels unique
    new_labels <- make.unique(base_labels, sep = "_")

    # apply the new labels to the dataframe
    tax_split[is_empty_or_na, t] <- new_labels

    # if genus we leave it blank when unknown and include family_sp. as species entry
    if (t == "Genus") {
        tax_split <- tax_split %>%
            mutate(
                empty_genus = if_else(
                    str_detect(Genus, "uncultured|unknown"),
                    "",
                    Genus
                )
            ) %>%
            mutate(
                fam_sp = if_else(
                    str_detect(Genus, "uncultured|unknown"),
                    paste(Family, "sp.", sep = "_"), Species
                )
            ) %>%
            mutate(Genus = empty_genus) %>%
            mutate(Species = fam_sp) %>%
            select(-c(fam_sp, empty_genus))
    }

    # if species is uncultured add genus in front as is already case with species annotated sequences
    if (t == "Species") {
        tax_split <- tax_split %>%
            mutate(
                genus_species = if_else(
                    str_detect(Species, "uncultured|metagenome|unidentified|unknown"),
                    paste(Genus, Species, sep = "_"),
                    Species
                )
            ) %>%
            mutate(Species = genus_species) %>%
            select(-genus_species)
    }
    # Make each "uncultured" annotation unique by adding an integer.
    tax_split[, t][grep("uncultured|metagenome|sp\\.", tax_split[, t])] <- make.unique(tax_split[, t][grep("uncultured|metagenome|sp\\.", tax_split[, t])], "_")
}

# remove spaces and *__ from names
tax_split <- tax_split %>%
    mutate_all(~ gsub(".__", "", .)) %>%
    mutate_all(~ gsub(" ", "", .))

# get counts table
otu_count <- otu_tax_edit_format %>%
    select(where(is.numeric))

# make it into a dataframe and name the rows with OTU and row number as with tax
otu_count <- as.data.frame(otu_count)
rownames(otu_count) <- paste0("OTU", seq_len(nrow(otu_count)))

# create phyloseq obj
otu <- phyloseq::otu_table(otu_count, taxa_are_rows = TRUE)
tax <- phyloseq::tax_table(as.matrix(tax_split))
sampledata <- phyloseq::sample_data(metadata)

physeq <- phyloseq(otu, tax, sampledata)

# filter physeq object
# get total reads
total_reads <- sum(unlist(as.data.frame(otu_table(physeq))))

# get the amount of reads corresponding to 0.1% of total reads
threshold_count0.1 <- total_reads * 0.1 / 100

# filter physeq object på threshold
ps <- prune_taxa(taxa_sums(physeq) > threshold_count0.1, physeq)

# split physeq into two based on sampling stage
physeq_filt_0.1_list <- split_physeq_in_two(
    physeq_obj = ps,
    feature = "sampling_stage",
    value1 = "L21",
    value2 = "G109"
)
ps_lactation <- physeq_filt_0.1_list[[1]]
ps_gestation <- physeq_filt_0.1_list[[2]]

# save the three phyloseq objects
save(ps, ps_gestation, ps_lactation,
    file = "data/physeq_objects_filt0.1.RData"
)
