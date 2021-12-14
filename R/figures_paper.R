source("R/load_data.R")
edna_leray_EPT <- 
  edna_leray_EPT %>% 
  filter(str_detect(site_code, "_U_"), count > 0)

edna_leese_EPT <- 
  edna_leese_EPT %>% 
  filter(str_detect(site_code, "_U_"), count > 0)

kicknet_EPT <-
  kicknet_EPT %>% 
  filter(str_detect(site_code, "_U_"), nr_ind > 0)


# Map
stream_location <- sites_meta %>%
  filter(site_location == "Upstream") %>%
  group_by(site_community) %>%
  summarise(x_coord = mean(x_coord),
            y_coord = mean(y_coord)) %>% 
  sf::st_as_sf(coords = c("x_coord", "y_coord"), crs = 21781)

sf::read_sf("data/Switzerland_shapefile/swissBOUNDARIES3D_1_3_TLM_LANDESGEBIET.shp") %>% 
  ggplot() +
  geom_sf(color = NA, fill = "grey") +
  geom_sf(data = stream_location) +
  theme_bw() +
  theme(legend.position = "none") +
  ggsflabel::geom_sf_text_repel(aes(label = site_community), color = "black",
                                data = stream_location) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.15) +
  ggspatial::annotation_north_arrow(location = "tl",   width = unit(1, "cm"))




# Figure 1
# Need object ref_G_EPT from AN_checks.R + reference database
list(CH_list = unique(EPT_list_CH$CH_list_species) %>% na.omit(),
     Kicknet = unique(kicknet_EPT$CH_list_species) %>% na.omit(),
     Ref_G = unique(ref_G_EPT$CH_list_species) %>% na.omit(),
     eDNA_leray = unique(edna_leray_EPT$CH_list_species) %>% na.omit(),
     eDNA_leese = unique(edna_leese_EPT$CH_list_species) %>% na.omit()) %>%
  UpSetR::fromList() %>% 
  rename(`EPT Switzerland` = CH_list,
         `Reference database` = Ref_G,
         `mlCOIintF/HCO2198` = eDNA_leray,
         `fwhF2/EPTDr2n` = eDNA_leese) %>% 
  UpSetR::upset(order.by = "freq",
                sets.x.label = "Number of species",
                mainbar.y.label = "Intersection (Number of species)",
                text.scale = 1.5)

#Simplified version
# list(CH_list = unique(EPT_list_CH$CH_list_species) %>% na.omit(),
#      Kicknet = unique(kicknet_EPT$CH_list_species) %>% na.omit(),
#      Ref_G = unique(ref_G_EPT$CH_list_species) %>% na.omit(),
#      eDNA_leray = unique(edna_leray_EPT$CH_list_species) %>% na.omit(),
#      eDNA_leese = unique(edna_leese_EPT$CH_list_species) %>% na.omit()) %>%
#   UpSetR::fromList() %>% 
#   filter(Kicknet > 0 | eDNA_leray > 0 | eDNA_leese > 0) %>% 
#   rename(`EPT Switzerland` = CH_list,
#          `Reference database` = Ref_G,
#          `eDNA (Leray)` = eDNA_leray,
#          `eDNA (Leese)` = eDNA_leese) %>% 
#   UpSetR::upset(order.by = "freq", text.scale = 2)


# Figure 2

edna_leese_kicknet_EPTs_frac <- 
  edna_leese_EPT %>% 
  left_join(sites_meta) %>%
  group_by(site_community, species_coal) %>% 
  summarise(count = sum(count)) %>% 
  full_join(group_by(left_join(kicknet_EPT, sites_meta), site_community, CH_list_species) %>% summarise(nr_ind = sum(nr_ind)),
            by = c("site_community" = "site_community", "species_coal" = "CH_list_species")) %>% 
  replace_na(replace = list(count = 0, nr_ind = 0)) %>%
  mutate(fraction = case_when(
    count > 0 & nr_ind > 0 ~ "Both",
    count > 0 & nr_ind == 0 ~ "eDNA",
    count == 0 & nr_ind > 0 ~ "Kicknet",
    TRUE ~ "None",
  )) %>% 
  filter(fraction != "None")

edna_leray_kicknet_EPTs_frac <- 
  edna_leray_EPT %>% 
  left_join(sites_meta) %>% 
  group_by(site_community, species_coal) %>% 
  summarise(count = sum(count)) %>% 
  full_join(group_by(left_join(kicknet_EPT, sites_meta), site_community, CH_list_species) %>% summarise(nr_ind = sum(nr_ind)),
            by = c("site_community" = "site_community", "species_coal" = "CH_list_species")) %>% 
  replace_na(replace = list(count = 0, nr_ind = 0)) %>%
  mutate(fraction = case_when(
    count > 0 & nr_ind > 0 ~ "Both",
    count > 0 & nr_ind == 0 ~ "eDNA",
    count == 0 & nr_ind > 0 ~ "Kicknet",
    TRUE ~ "None",
  )) %>% 
  filter(fraction != "None")


bind_rows(edna_leray_kicknet_EPTs_frac,
          edna_leese_kicknet_EPTs_frac,
          .id = "method") %>% 
  group_by(method, site_community, fraction) %>% 
  count() %>% 
  mutate(n = ifelse(method == "1", n * -1, n)) %>% 
  ggplot() +
  geom_col(aes(forcats::fct_rev(site_community), n, fill = fraction)) +
  geom_point(aes(a, z), data = tibble(a = "Aadorf", z = -53, method = "1"), color = "white", alpha = 0) +
  ggpol::facet_share(~method, dir = "h", reverse_num = FALSE, scales = "free",
                     labeller = as_labeller(c("1" = "mlCOIintF/HCO2198", "2" = "fwhF2/EPTDr2n"))) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "bottom", axis.title.y = element_blank()) +
  ylab("Number of taxa") +
  scale_y_continuous(breaks = c(-40, -20, 0, 20, 40),
                     labels = c(40, 20, 0, 20, 40)) +
  scale_fill_brewer(type = "qual", palette = 7)


bind_rows(edna_leray_kicknet_EPTs_frac,
          edna_leese_kicknet_EPTs_frac,
          .id = "method") %>% 
  group_by(method, species_coal, fraction) %>% 
  count() %>% 
  group_by(method, species_coal) %>% 
  mutate(tot = sum(n)) %>%
  group_by(species_coal) %>%
  mutate(max_tot = max(tot)) %>% 
  filter(max_tot > 1) %>% 
  ggplot() +
  geom_col(aes(forcats::fct_rev(species_coal), n, fill = fraction)) +
  facet_wrap(vars(method), labeller = as_labeller(c("1" = "mlCOIintF/HCO2198", "2" = "fwhF2/EPTDr2n"))) +
  coord_flip() +
  ylab("Number of sites") +
  theme(axis.title.y = element_blank()) +
  scale_fill_brewer(type = "qual", palette = 7)


alpha_richness_leray <- full_join(
  edna_leray_kicknet_EPTs_frac %>% 
    group_by(site_community) %>% 
    count(fraction) %>% 
    filter(fraction != "eDNA") %>% 
    summarise(KN = sum(n)),
  
  edna_leray_kicknet_EPTs_frac %>% 
    group_by(site_community) %>% 
    count(fraction) %>% 
    filter(fraction != "Kicknet") %>% 
    summarise(eDNA = sum(n))
) %>% 
  replace_na(list(KN = 0, eDNA = 0))

alpha_richness_leese <- full_join(
  edna_leese_kicknet_EPTs_frac %>% 
    group_by(site_community) %>% 
    count(fraction) %>% 
    filter(fraction != "eDNA") %>% 
    summarise(KN = sum(n)),
  edna_leese_kicknet_EPTs_frac %>% 
    group_by(site_community) %>% 
    count(fraction) %>% 
    filter(fraction != "Kicknet") %>% 
    summarise(eDNA = sum(n))
) %>% 
  replace_na(list(KN = 0, eDNA = 0))

t.test(alpha_richness_leray$KN, alpha_richness_leray$eDNA)
t.test(alpha_richness_leese$KN, alpha_richness_leese$eDNA)
t.test(alpha_richness_leray$eDNA, alpha_richness_leese$eDNA, paired = TRUE)
sd(alpha_richness_leray$eDNA)
sd(alpha_richness_leese$eDNA)


# Figure 3
edna_leray_EPT_richness <- edna_leray_EPT %>% 
  group_by(site_code, species_coal) %>% 
  summarise(count = sum(count)) %>% 
  spread_cdm(site_code, species_coal, count, fill.missing = 0) %>% 
  vegan::specnumber() %>% 
  enframe("site_code", "Leray")

edna_leese_EPT_richness <- edna_leese_EPT %>% 
  group_by(site_code, species_coal) %>% 
  summarise(count = sum(count)) %>% 
  spread_cdm(site_code, species_coal, count, fill.missing = 0) %>% 
  vegan::specnumber() %>% 
  enframe("site_code", "Leese")

kicknet_EPT_richness <- kicknet_EPT %>% 
  group_by(site_code, CH_list_species) %>% 
  summarise(count = sum(nr_ind)) %>% 
  spread_cdm(site_code, CH_list_species, count, fill.missing = 0) %>% 
  vegan::specnumber() %>% 
  enframe("site_code", "KN")


EPT_richness <- 
  kicknet_EPT_richness %>% 
  left_join(edna_leray_EPT_richness, by = "site_code") %>% 
  left_join(edna_leese_EPT_richness, by = "site_code") %>% 
  left_join(sites_meta, by = "site_code") %>% 
  filter(site_location == "Upstream") %>% 
  replace_na(list(Leray = 0, Leese = 0)) %>% 
  group_by(site_community, site_location) %>% 
  summarise(KN = mean(KN),
            Leray = mean(Leray),
            Leese = mean(Leese)) %>% 
  select(-site_location) %>% 
  left_join(model_kaelin, by = "site_community") %>% 
  select(-Step_IBCH, -Lasso_IBCH, -Step_EPT) %>%
  ungroup()


EPT_richness %>% 
  select(-site_community) %>% 
  set_colnames(c("Kicknet", "mlCOIintF/HCO2198", "fwhF2/EPTDr2n", "Lasso model")) %>% 
  GGally::ggpairs(
    lower = list(continuous = GGally::wrap(lowerFn, method = "lm")),
    upper = list(continuous = GGally::wrap("cor", size = 6))
  ) +
  theme_light()


mat_test <- EPT_richness %>% 
  select(-site_community)

cor.test(mat_test$Lasso_EPT, mat_test$Leese)
cor.test(mat_test$Leese, mat_test$KN)


# SI

# Rarefaction curves Leese
edna_leese_EPT %>% 
  left_join(sites_meta) %>% 
  group_by(site_community, species_coal) %>% 
  summarise(count = sum(count)) %>% 
  spread_cdm(site_community, species_coal, count, fill.missing = 0) -> zzz

ggg <- vegan::rarecurve(zzz) %>% 
  tidy_rarecurve(rownames(zzz))

si_1 <- ggplot(ggg) +
  geom_line(aes(samp_size, species, group = site, color = site)) +
  xlab("Sample size (reads)") +
  ylab("Species")

# Rarefaction curves Leray
edna_leray_EPT %>% 
  left_join(sites_meta) %>% 
  group_by(site_community, species_coal) %>% 
  summarise(count = sum(count)) %>% 
  spread_cdm(site_community, species_coal, count, fill.missing = 0) -> zzz

ggg <- vegan::rarecurve(zzz) %>% 
  tidy_rarecurve(rownames(zzz))

si_2 <- ggplot(ggg) +
  geom_line(aes(samp_size, species, group = site, color = site)) +
  xlab("Sample size (reads)") +
  ylab("Species")


# Rarefaction curves Kicknet
kicknet_EPT %>% 
  left_join(sites_meta) %>% 
  group_by(site_community, CH_list_species) %>% 
  summarise(count = sum(nr_ind)) %>% 
  spread_cdm(site_community, CH_list_species, count, fill.missing = 0) -> zzz

ggg <- vegan::rarecurve(zzz) %>% 
  tidy_rarecurve(rownames(zzz))

si_3 <- ggplot(ggg) +
  geom_line(aes(samp_size, species, group = site, color = site)) +
  xlab("Sample size (individuals)") +
  ylab("Species")



edna_both_EPT_richness <- edna_leese_EPT %>% 
  bind_rows(edna_leray_EPT) %>% 
  group_by(site_code, species_coal) %>% 
  summarise(count = sum(count)) %>% 
  spread_cdm(site_code, species_coal, count, fill.missing = 0) %>% 
  vegan::specnumber() %>% 
  enframe("site_code", "Merged eDNA")


si_4 <- kicknet_EPT_richness %>% 
  left_join(edna_both_EPT_richness, by = "site_code") %>% 
  left_join(sites_meta, by = "site_code") %>% 
  filter(site_location == "Upstream") %>% 
  replace_na(list(`Merged eDNA` = 0)) %>% 
  group_by(site_community, site_location) %>% 
  summarise(KN = mean(KN),
            `Merged eDNA` = mean(`Merged eDNA`)) %>% 
  select(-site_location) %>% 
  left_join(model_kaelin, by = "site_community") %>% 
  select(-Step_IBCH, -Lasso_IBCH, -Step_EPT) %>%
  ungroup() %>% 
  select(-site_community) %>% 
  set_colnames(c("Kicknet", "Merged eDNA", "Lasso model")) %>% 
  GGally::ggpairs(
    lower = list(continuous = GGally::wrap(lowerFn, method = "lm")),
    upper = list(continuous = GGally::wrap("cor", size = 6))
  ) +
  theme_light()


# Compile Supplementary Material
rmarkdown::render("Rmd/supplementary_info.Rmd", output_dir = "results", output_file = "supplementary_info.pdf")
browseURL("results/supplementary_info.pdf")
