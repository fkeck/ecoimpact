
#Unzip data/Ref_mothur_reftax_val.zip
#And edit the file path accordingly
ref_G <- read_delim("/home/ecoadmin/Documents/postdoc/ReferenceG_v290420.fa/ReferenceG_mothur_reftax_val",
                        delim = ";",
                        col_names = c("ID_Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Xtra")) %>% 
  separate(ID_Domain, into = c("ID", "Domain"), sep = " ") %>% 
  select(-Xtra) %>% 
  mutate(Species = str_replace_all(Species, "_", " "))
  

# EPT
ref_G_EPT <- ref_G %>% 
  filter(Order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera")) %>% 
  left_join(select(EPT_list_CH, species, CH_list_species), by = c("Species" = "species"))

tibble(Species = unique(c(EPT_list_CH$CH_list_species,
                          kicknet_EPT$CH_list_species,
                          edna_leray_EPT$species_coal,
                          edna_leese_EPT$species_coal,
                          ref_G_EPT$CH_list_species))) %>% 
  mutate(CH_list = Species %in% EPT_list_CH$CH_list_species,
         Kicknet = Species %in% kicknet_EPT$CH_list_species,
         Ref_G = Species %in% ref_G_EPT$CH_list_species,
         eDNA_leray = Species %in% edna_leray_EPT$CH_list_species,
         eDNA_leese = Species %in% edna_leese_EPT$CH_list_species) %>% 
  write_csv("results/species_bool_datasets.csv")


list(CH_list = unique(EPT_list_CH$CH_list_species) %>% na.omit(),
     Kicknet = unique(kicknet_EPT$CH_list_species) %>% na.omit(),
     Ref_G = unique(ref_G_EPT$CH_list_species) %>% na.omit(),
     eDNA_leray = unique(edna_leray_EPT$CH_list_species) %>% na.omit(),
     eDNA_leese = unique(edna_leese_EPT$CH_list_species) %>% na.omit()) %>%
  UpSetR::fromList() %>% 
  UpSetR::upset(order.by = "freq")



# IBCH
ref_G_IBCH <- ref_G %>% 
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


tibble(Family = unique(c(IBCH_list_CH$CH_list,
                         kicknet_IBCH$family,
                         edna_leray_IBCH$Family,
                         edna_leese_IBCH$Family,
                         ref_G_IBCH$Family))) %>% 
  mutate(Kicknet = Family %in% kicknet_IBCH$family,
         Ref_G = Family %in% ref_G_IBCH$Family,
         eDNA_leray = Family %in% edna_leray_IBCH$Family,
         eDNA_leese = Family %in% edna_leese_IBCH$Family) %>% 
  write_csv("results/families_bool_datasets.csv")


list(CH_list = unique(IBCH_list_CH$CH_list) %>% na.omit(),
     Kicknet = unique(kicknet_IBCH$family) %>% na.omit(),
     Ref_G = unique(ref_G_IBCH$Family) %>% na.omit(),
     eDNA_leray = unique(edna_leray_IBCH$Family) %>% na.omit(),
     eDNA_leese = unique(edna_leese_IBCH$Family) %>% na.omit()) %>%
  UpSetR::fromList() %>% 
  UpSetR::upset(order.by = "freq")



ref_seq <- read_fasta("/home/ecoadmin/Documents/postdoc/ReferenceG_v290420.fa/test_mothur/ReferenceG_mothur.fasta")

ref_G_EPT %>% 
  left_join(as_tibble.bioseq_dna(ref_seq, label = "ID")) %>% 
  filter(!is.na(sequence),
         CH_list_species %in% unique(EPT_list_CH$CH_list_species) %>% na.omit()) %>% 
  mutate(
    Leese = seq_crop_pattern(sequence,
                             pattern_in = NULL,
                             pattern_out = seq_reverse(seq_complement(dna("CAAACAAATARDGGTATTCGDTY"))),
                             include_patterns = F)
  )


