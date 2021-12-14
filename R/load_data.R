library(tidyverse)
library(tidymodels)
library(bioseq)
library(flexitarian)
library(magrittr)

theme_set(theme_classic())
plot_list <- list()
stats_list <- list()

sites_meta <- read_csv("data/sites_meta.csv")



landuse <- read_csv("data/site_landuse.csv") %>% 
  select(starts_with(c("site_community", "fraction_")))

landuse_frac <- c("fraction_arable",
                  "fraction_orchard",
                  "fraction_pasture",
                  "fraction_forested",
                  "fraction_urban",
                  "fraction_agri_intense")

EPT_list_CH <- read_csv("data/EPT_list_CH.csv")
IBCH_list_CH <- read_csv("data/IBCH_list_CH.csv")

kicknet <- read_csv("data/kicknet_data.csv")

kicknet_EPT <- kicknet %>% 
  filter(identification == "EPT")

kicknet_IBCH <- kicknet %>% 
  filter(identification == "IBCH",
         family %in% IBCH_list_CH$CH_list)

model_kaelin <- read_csv("data/model_pred_Kaelin.csv")

edna_leray_files <- read_csv("data/leray_files.csv")

edna_leray_id <- read_delim("data/ASV_to_id.ReferenceG_mothur_ref.wang.taxonomy",
           delim = ";",
           col_names = c("ID_Domain", "Phylum", "Class",
                         "Order", "Family", "Genus",
                         "Species", "Xtra")) %>% 
  separate(ID_Domain, into = c("ID", "Domain"), sep = "\t") %>% 
  select(-Xtra) %>% 
  mutate(across(Domain:Species, str_replace_all, "_", " ")) %>% 
  mutate(across(everything(), str_remove_all, "\\(.+\\)")) %>% 
  left_join(read_fasta("data/ASV_to_id.fasta") %>%
              as_tibble.bioseq_dna(label = "ID", sequence = "DNA_seq"))


edna_leray <- read_csv("data/leray_sequence_table_nochim.csv") %>% 
  rename(DNA_seq = DNA_SEQ) %>% 
  mutate(DNA_seq = as_dna(DNA_seq))


# Filter out ASV with codon stop
edna_leray <- edna_leray %>% 
  mutate(AA_seq = seq_translate(DNA_seq, code = 5, codon_frame = 2)) %>% 
  filter(!seq_detect_pattern(AA_seq, "\\*|X")) %>% 
  select(-AA_seq) %>% 
  pivot_longer(!DNA_seq, names_to = "seq_root_file", values_to = "count") %>% 
  left_join(edna_leray_files %>% distinct(seq_root_file, site_code)) %>% 
  select(DNA_seq, site_code, count)

# Filter out negative control
edna_leray <- edna_leray %>%
  filter(str_starts(site_code, "NEG_")) %>% 
  left_join(edna_leray_id) %>% 
  mutate(prop = count/sum(count)) %>% 
  filter(prop > 0.001) %>% 
  anti_join(edna_leray, ., by = "DNA_seq") %>% 
  filter(!str_starts(site_code, "NEG_"))


# Filter out effluent samples
edna_leray <- edna_leray %>%
  filter(!str_detect(site_code, "_E$"))

# Flush empty ASV
edna_leray <- edna_leray %>% 
  group_by(DNA_seq) %>% 
  mutate(tot_count = sum(count)) %>% 
  filter(tot_count > 0) %>% 
  select(-tot_count)

edna_leray_id <- edna_leray_id %>% 
  filter(DNA_seq %in% unique(edna_leray$DNA_seq))

# Filter ASV with < 10 occurrences in the whole dataset
# edna_leray <- edna_leray %>%
#   group_by(DNA_seq) %>% 
#   mutate(total_count = sum(count)) %>% 
#   filter(total_count > 10) %>% 
#   select(-total_count)


# Filter ASV unclassified
# edna_leray <- edna_leray %>% 
#   left_join(edna_leray_id) %>% 
#   filter(Class != "Eukaryota unclassified")


# Normalize using proportions
# edna_leray <- edna_leray %>% 
#   group_by(site_code) %>% 
#   mutate(prop = count/sum(count),
#          pa = ifelse(count > 0, 1L, 0L),
#          .after = count)


# EPT
edna_leray_EPT <- edna_leray %>% 
  left_join(edna_leray_id) %>%
  filter(Order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera")) %>%
  left_join(EPT_list_CH, by = c("Species" = "species")) %>% 
  mutate(genus_coal = coalesce(CH_list_genus, Genus),
         species_coal = coalesce(CH_list_species, Species), .after = ID)

# IBCH
edna_leray_IBCH <- edna_leray %>% 
  left_join(edna_leray_id) %>%
  mutate(Family = ifelse(Family %in% c("Limoniidae", "Pedicidae"), "Limoniidae_Pedicidae", Family)) %>% 
  mutate(Family = ifelse(Family %in% c("Anthomyiidae", "Muscidae"), "Anthomyiidae_Muscidae", Family)) %>% 
  mutate(Family = ifelse(Phylum == "Porifera", "PORIFERA", Family)) %>% 
  mutate(Family = ifelse(Phylum == "Cnidaria", "CNIDARIA", Family)) %>% 
  mutate(Family = ifelse(Phylum == "Bryozoa", "BRYOZOA", Family)) %>% 
  mutate(Family = ifelse(Phylum == "Nematoda", "NEMATHELMINTHES", Family)) %>% 
  mutate(Family = ifelse(Order %in% c("Enchytraeida", "Haplotaxida", "Lumbriculida", "Moniligastrida"), "Oligochaeta", Family)) %>%
  mutate(Family = ifelse(Class == "Branchiopoda", "Branchiopoda", Family)) %>%
  mutate(Family = ifelse(Order == "Hymenoptera", "Hymenoptera", Family)) %>%
  mutate(Family = ifelse(Order == "Lepidoptera", "Lepidoptera", Family)) %>%
  filter(Family %in% IBCH_list_CH$CH_list)











edna_leese_files <- read_csv("data/leese_files.csv")

edna_leese_id <- read_delim("data/ASV_to_id_Leese_RefG.wang.taxonomy",
                            delim = ";",
                            col_names = c("ID_Domain", "Phylum", "Class",
                                          "Order", "Family", "Genus",
                                          "Species", "Xtra")) %>% 
  separate(ID_Domain, into = c("ID", "Domain"), sep = "\t") %>% 
  select(-Xtra) %>% 
  mutate(across(Domain:Species, str_replace_all, "_", " ")) %>% 
  mutate(across(everything(), str_remove_all, "\\(.+\\)")) %>% 
  left_join(read_fasta("data/ASV_to_id_Leese.fasta") %>%
              as_tibble.bioseq_dna(label = "ID", sequence = "DNA_seq"))


edna_leese <- read_csv("data/leese_sequence_table_nochim.csv") %>% 
  rename(DNA_seq = DNA_SEQ) %>% 
  mutate(DNA_seq = as_dna(DNA_seq))


# Filter out ASV with codon stop
edna_leese <- edna_leese %>% 
  mutate(AA_seq = seq_translate(DNA_seq, code = 5, codon_frame = 2)) %>% 
  filter(!seq_detect_pattern(AA_seq, "\\*|X")) %>% 
  select(-AA_seq) %>% 
  pivot_longer(!DNA_seq, names_to = "seq_root_file", values_to = "count") %>% 
  left_join(edna_leese_files %>% distinct(seq_root_file, site_code)) %>% 
  select(DNA_seq, site_code, count)

# Filter out negative control
edna_leese <- edna_leese %>%
  filter(str_starts(site_code, "PCR_NEG_")) %>% 
  left_join(edna_leese_id) %>% 
  mutate(prop = count/sum(count)) %>% 
  filter(prop > 0.001) %>% 
  anti_join(edna_leese, ., by = "DNA_seq") %>% 
  filter(!str_starts(site_code, "PCR_NEG_"))

# Filter out bad quality samples
 edna_leese <- edna_leese %>%
   filter(!site_code %in% c("BUT_U_2"))

# Filter out effluent samples
edna_leese <- edna_leese %>%
  filter(!str_detect(site_code, "_E$"))


# FILTER EX, FILTERS AND POSITIVES 
edna_leese <- edna_leese %>%
  filter(!str_starts(site_code, "EX_|FILTER_|PCR_POS|UND"))


# Flush empty ASV
edna_leese <- edna_leese %>% 
  group_by(DNA_seq) %>% 
  mutate(tot_count = sum(count)) %>% 
  filter(tot_count > 0) %>% 
  select(-tot_count)

edna_leese_id <- edna_leese_id %>% 
  filter(DNA_seq %in% unique(edna_leese$DNA_seq))


# EPT
edna_leese_EPT <- edna_leese %>% 
  left_join(edna_leese_id) %>%
  filter(Order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera")) %>%
  left_join(EPT_list_CH, by = c("Species" = "species")) %>% 
  mutate(genus_coal = coalesce(CH_list_genus, Genus),
         species_coal = coalesce(CH_list_species, Species), .after = ID)

# IBCH
edna_leese_IBCH <- edna_leese %>% 
  left_join(edna_leese_id) %>%
  mutate(Family = ifelse(Family %in% c("Limoniidae", "Pedicidae"), "Limoniidae_Pedicidae", Family)) %>% 
  mutate(Family = ifelse(Family %in% c("Anthomyiidae", "Muscidae"), "Anthomyiidae_Muscidae", Family)) %>% 
  mutate(Family = ifelse(Phylum == "Porifera", "PORIFERA", Family)) %>% 
  mutate(Family = ifelse(Phylum == "Cnidaria", "CNIDARIA", Family)) %>% 
  mutate(Family = ifelse(Phylum == "Bryozoa", "BRYOZOA", Family)) %>% 
  mutate(Family = ifelse(Phylum == "Nematoda", "NEMATHELMINTHES", Family)) %>% 
  mutate(Family = ifelse(Order %in% c("Enchytraeida", "Haplotaxida", "Lumbriculida", "Moniligastrida"), "Oligochaeta", Family)) %>%
  mutate(Family = ifelse(Class == "Branchiopoda", "Branchiopoda", Family)) %>%
  mutate(Family = ifelse(Order == "Hymenoptera", "Hymenoptera", Family)) %>%
  mutate(Family = ifelse(Order == "Lepidoptera", "Lepidoptera", Family)) %>%
  filter(Family %in% IBCH_list_CH$CH_list)



lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "blue") +
    geom_smooth(method = method, color = "red", ...)
  p
}




#' Convert vegan rarecurve results to tidy format
#'
#' @param x a list of rarefy results returned by rarecurve
#' @param sites a vector of site names
#'
#' @return A data.frame in "long format".
#'
tidy_rarecurve <- function(x, sites) {
  long_list <- mapply(function(x, sites) {
    data.frame(
      site = sites,
      samp_names = names(x),
      samp_size = attr(x, "Subsample"),
      species = x)},
    x = x, sites = sites, SIMPLIFY = FALSE)
  
  do.call(rbind, long_list)
}
