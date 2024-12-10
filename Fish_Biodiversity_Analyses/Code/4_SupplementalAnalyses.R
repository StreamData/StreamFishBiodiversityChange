library(geofacet)

##Supplemental analyses:

#test thermal preference####
#compare Ctmax from Bayat et al. database to thermal preference from fishlife

ctmax_fish = read.csv("./Data/heattol_fish_ctmax.csv")

unique(ctmax_fish$endpoint)

ctmax_fishsub = ctmax_fish %>%
  filter(metric == "ctmax") %>%
  mutate(species = gsub(" ",".",species)) %>%
  dplyr::select(species, metric, tmax, n, error:origin,continent_f ) %>%
  filter(species %in% row.names(FishLifeTraits)) %>%
  filter(!is.na(tmax))

tempsub = FishLifeTraits %>%
  mutate(species = row.names(FishLifeTraits)) %>%
  dplyr::select(species, temperature) %>%
  filter(species %in% ctmax_fishsub$species)

ctmax_wtemps = ctmax_fishsub %>%
  left_join(tempsub) %>%
  mutate(acclim_temp = as.numeric(acclim_temp)) 

mctmax = lmer(tmax ~ temperature + acclim_temp + (acclim_temp|species),
     data = ctmax_wtemps)

car::Anova(mctmax, type = 2, test = "F")
performance::r2(mctmax)

ctmaxfigdat <- data.frame(ggeffects::ggemmeans(mctmax, terms = c("temperature[all]")))

ggplot(ctmaxfigdat, aes(x = x, y = predicted, group = group))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),alpha = 0.25,color = NA)+
  geom_line(size = 1)+
  ylab(expression("CT"[max] ~(degree*C*", Bayat et al.")))+
  xlab(expression("Thermal perference " (degree*C*", FishLife")))+
  theme_bw()+
  theme(
    axis.text = element_text(color = "black",
                             size = 9),
    axis.title = element_text(size = 11))

ggsave("./Figures/Ctmax.jpg", dpi = 300, width = 4, height = 3.5, unit = "in")

#Figure S3----
##Non-rarefied Species Richness####
#use glmmTMB, because glmer.nb would not converge
fishdat$SppRichness = rowSums(fishdat[,24:412])

m_rawrich <- glmmTMB(SppRichness ~ Year*HUC2+Agency + StreamOrder + SampleTypeCode +
                   log(PredictedWettedWidth_km) + poly(log(WholeConductivity),2,raw = TRUE) +
                   Year*wt_pred_new+
                   Landuse+
                   (1|SiteNumber) + (1|YearC),
                 data = fishdat,
                 family = "poisson",
                 control=glmmTMBControl(optimizer=optim,
                                        optArgs=list(method="BFGS")))

car::Anova(m_rawrich, type = 3)
performance::r2(m_rawrich)


frrdat <- data.frame(ggeffects::ggemmeans(m_rawrich, 
                                          terms = c("Year[all]","wt_pred_new[13.0,20.5,25.6]"),
                                          weights = 'proportional'))
frrdat$group = as.factor(frrdat$group )

test(emtrends(m_rawrich, ~wt_pred_new, var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = 'proportional'))
# wt_pred_new Year.trend      SE  df z.ratio p.value
# 13.0  -1.32e-02 0.00610 Inf  -2.158  0.0309
# 20.5   1.59e-05 0.00155 Inf   0.010  0.9918
# 25.6   8.98e-03 0.00369 Inf   2.431  0.0151

rawrf1 = ggplot(frrdat, aes(x = x, y = predicted, group = group))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_color_manual(values = rev(c("red","violet","blue")),
                     labels = c("Cold","Intermediate","Warm"),
                     name = "Past temperature regime")+
  scale_fill_manual(values = rev(c("red","violet","blue")),
                    labels = c("Cold","Intermediate","Warm"),
                    name = "Past temperature regime")+
  scale_y_continuous(limits = c(3.75,16),
                     breaks = c(5,10,15))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  ylab("Non-rarefied species richness")+
  xlab("Year")+
  theme_bw()+
  theme(
        legend.position = "bottom",
        legend.margin = margin(-5,10,-5,-10),
        axis.text = element_text(color = "black",
                                 size = 9),
        axis.title = element_text(size = 11))

##Evenness####

S_pie = diversity(fishdatCPUE[,24:405], index = "invsimpson")

##sites with 0 species are given NA values
fishdatCPUE$Spie = ifelse(is.infinite(S_pie),
                           NA,
                           S_pie)

m_Spie <- lmer(log(Spie) ~ Year*HUC2+Agency +
                       StreamOrder + SampleTypeCode +
                       SANDCAT+
                       log(PredictedWettedWidth_m) +
                       poly(log(WholeConductivity),2,raw = TRUE) +
                       Year*wt_pred_new +
                       Landuse+
                       (1|SiteNumber) + (1|YearC),
                     data = fishdatCPUE)

car::Anova(m_Spie, type = 3)

performance::r2(m_Spie)

##Significant interaction; but significant evenness differences
##cold and intermediate streams significantly decreasing
test(emtrends(m_Spie, ~wt_pred_new, var = "Year", rg.limit = 1000000, lmerTest.limit = 10,
              pbkrtest.limit = 10,
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = 'proportional'))
# wt_pred_new Year.trend      SE  df z.ratio p.value
# 13.0  -0.009116 0.00364 Inf  -2.504  0.0123
# 20.5  -0.000489 0.00176 Inf  -0.277  0.7814
# 25.6   0.005378 0.00262 Inf   2.051  0.0402


fspiedat <- data.frame(ggeffects::ggemmeans(m_Spie,
                                            terms = c("Year[all]","wt_pred_new[13.0,20.5,25.6]"),
                                            weights = 'proportional'))
fspiedat$group = as.factor(fspiedat$group )

p1 = ggplot(fspiedat, aes(x = x, y = predicted, group = group))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_color_manual(values = rev(c("red","violet","blue")),
                     labels = c("Cold","Intermediate","Warm"),
                     name = "Past temperature regime")+
  scale_fill_manual(values = rev(c("red","violet","blue")),
                    labels = c("Cold","Intermediate","Warm"),
                    name = "Past temperature regime")+
  # scale_y_continuous(limits = c(1.5,6))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  ylab(expression(Evenness~(S[PIE])))+
  xlab("Year")+
  theme_bw()+
  theme(
    legend.position = "none",
    legend.margin = margin(-5,10,-5,-10),
    axis.text = element_text(color = "black",
                             size = 9),
    axis.title = element_text(size = 11))

legendrp = get_plot_component(rawrf1 + theme(legend.box.margin = margin(c(0,0,0,0))),
                              'guide-box-bottom', return_all = TRUE)

suppfig2 = plot_grid(plot_grid(rawrf1 + theme(legend.position = "none"), p1,
                    align = "hv",
                    nrow = 1,
                    axis = "l",
                    labels = c("A","B"),
                    label_x = .075),
          legendrp,
          ncol = 1,
          rel_heights = c(1,.1))

ggsave("./Figures/FigS2.jpg", suppfig2,dpi = 600, width = 5.75, height = 3, unit = "in")



#Spatial analyses####
#Note: version of spmodel used here is the devbranch version, as of 8/15/24 the
# emmeans extensions had not been incorporated into the main branch of spmodel
library(spmodel)

m1.sp <- splm(log(TotalCPUE100m + 0.01) ~ Year*HUC2 + Agency +
                StreamOrder  +
                SANDCAT+
                SampleTypeCode + log(MethodEffort)+ log(ReachLengthFished_100m)+
                log(PredictedWettedWidth_m) +
                poly(log(WholeConductivity), 2, raw = TRUE) +
                Year*wt_pred_new +
                Year*Landuse,
              data = fishdatCPUE,
              random = ~ SiteNumber + YearC,
              spcov_type = "exponential",
              xcoord = Longitude_dd, ycoord = Latitude_dd
)

test(emtrends(m1.sp,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE  df z.ratio p.value
# 12.9   -0.01948 0.00997 Inf  -1.954  0.0507
# 19.8    0.00294 0.00567 Inf   0.519  0.6035
# 24.9    0.01952 0.00669 Inf   2.919  0.0035

m_rare.sp <- splm(log(RarefiedRichness) ~  Year*HUC2+Agency +
                    StreamOrder + SampleTypeCode +
                    log(PredictedWettedWidth_m) +
                    poly(log(WholeConductivity),2,raw = TRUE) +
                    Year*wt_pred_new+Landuse,
                  data = fishdat_rarefiedRich,
                  random = ~ SiteNumber + YearC,
                  spcov_type = "exponential",
                  xcoord = Longitude_dd, ycoord = Latitude_dd
)

anova(m_rare.sp)



test(emtrends(m_rare.sp,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE  df z.ratio p.value
# 12.9  -0.016525 0.00296 Inf  -5.579  <.0001
# 19.8  -0.006525 0.00159 Inf  -4.101  <.0001
# 24.9   0.000867 0.00207 Inf   0.419  0.6752


m1f.sp <- splm(SES ~  Year*HUC2 + Agency + StreamOrder + SampleTypeCode +
                 SANDCAT+
                 log(PredictedWettedWidth_km)+
                 poly(log(WholeConductivity), 2, raw = TRUE)+
                 wt_pred_new*Year+
                 Landuse,
               data = fishdat,
               random = ~ SiteNumber + YearC,
               spcov_type = "exponential",
               xcoord = Longitude_dd, ycoord = Latitude_dd
)

test(emtrends(m1f.sp,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE  df z.ratio p.value
# 12.9    0.01524 0.00758 Inf   2.011  0.0444
# 19.8    0.00197 0.00374 Inf   0.527  0.5983
# 24.9   -0.00784 0.00509 Inf  -1.540  0.1235

mLCBD.sp <- splm(LCBD_s2 ~  Year*HUC2 + Agency + StreamOrder + SampleTypeCode +
                   SANDCAT+
                   log(PredictedWettedWidth_km)+
                   poly(log(WholeConductivity), 2, raw = TRUE)+
                   wt_pred_new*Year+
                   Landuse,
                 data = fishdat_LCBD,
                 random = ~ SiteNumber + YearC,
                 spcov_type = "exponential",
                 xcoord = Longitude_dd, ycoord = Latitude_dd
)

test(emtrends(mLCBD.sp,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend       SE  df z.ratio p.value
# 12.9   2.30e-03 0.000817 Inf   2.812  0.0049
# 19.8   4.18e-05 0.000458 Inf   0.091  0.9273
# 24.9  -1.63e-03 0.000583 Inf  -2.790  0.0053


#life history supplemental analyses----

###Native, non-game compared to Non-native, game ####
OPE_sub = OPE_sub %>%
  mutate(species = ifelse(species == "Lethenteron appendix",
                          "Lampetra appendix",
                          ifelse(species == "Ericymba amplamala",
                                 "Notropis amplamala",
                                 ifelse(species == "Etheostoma chlorosomum",
                                        "Etheostoma chlorosoma",
                                        species))))

OPE_sublink = OPE_sub %>%
  mutate(species = gsub(" ", ".",species),
         genus = sub("\\..*", "", species)) %>%
  filter(genus %in% c("Cottus", "Moxostoma","Hybopsis", "Notropis")) %>%
  group_by(genus) %>%
  summarize(across(Equilibrium:Opportunist, ~mean(.)))  %>%
  mutate(freq = ifelse(genus %in% c("Hybopsis","Moxostoma"),
                       2,
                       ifelse(genus == "Notropis",
                              3,
                              1))) %>%
  uncount(freq)

OPE_sublink$species = c("Cottus chattahoochee",
                        "Hybopsis winchelli",
                        "Hybopsis rubrifrons",
                        "Moxostoma duquesnei",
                        "Moxostoma collapsum",
                        "Notropis percobromus",
                        "Notropis candidus",
                        "Notropis alborus")

OPE_subfin = bind_rows(OPE_sub,
                       OPE_sublink %>%
                         dplyr::select(-genus))


fish_nat_stat <- read.csv("./data/NonnativeHUC2.csv",
                          colClasses = list("HUC2" = "character"))
gamey = read.csv("./Data/GameFishDesignation.csv")


fish_nat_stat = fish_nat_stat %>%
  group_by(Scientific.Name) %>%
  slice(1) %>%
  dplyr::select(-HUC2)


gamey = gamey %>%
  filter(Species %in% OPE_subfin$species) %>%
  mutate(Species = gsub(" ",".",Species))

fishstats = OPE_subfin %>%
  mutate(species = gsub(" ",".",species)) %>%
  left_join(fish_nat_stat, join_by(species == "Scientific.Name")) %>%
  mutate(Native = ifelse(is.na(Native),
                         "Native",
                         Native)) %>%
  left_join(gamey, 
            join_by(species == Species)) %>%
  mutate(Status = ifelse(is.na(Status),
                         "Not Game",
                         Status)) %>%
  # mutate(NNG = ifelse(Status == "Game fish",
  #                     "Game fish",
  #                     "Not game"))
  
  mutate(NNG = ifelse(Status == "Game fish" | Native != "Native",
                      "Non-native, game",
                      "Native, non-game")) %>%
  mutate(across(Equilibrium:Opportunist, ~((. * (389 - 1)) + 0.5) / 389))

mO = glmmTMB(Opportunist ~ Native * Status,
             family = beta_family(link = "logit"),
             data = fishstats,
             control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
car::Anova(mO, type = 3)

mE = glmmTMB(Equilibrium ~ Native * Status,
             family = beta_family(link = "logit"),
             data = fishstats,
             control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
car::Anova(mE, type = 3)

mP = glmmTMB(Periodic ~ Native * Status,
             family = beta_family(link = "logit"),
             data = fishstats,
             control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
car::Anova(mP, type = 3)




#family-HUC2 trends in occupancy####

HUC2 <- read_sf("./Data/HUC2_ClipUS.shp")
HUC2 <- st_transform(HUC2,crs = 5070)

US_Bounds <- read_sf("./Data/UnitedStates.shp")
US_Bounds <- st_transform(US_Bounds,crs = 5070)

HUC2$Region = ifelse(HUC2$huc2 %in% c("02","03","08","12"),
                     "Coastal plains",
                     ifelse(HUC2$huc2 %in% c("04","05","07","09"),
                            "Midwest",
                            ifelse(HUC2$huc2 %in% c("13","14","15","16"),
                                   "Xeric",
                                   ifelse(HUC2$huc2 %in% c("10","11"),
                                          "Great plains",
                                          ifelse(HUC2$huc2 %in% c("17","18"),
                                                 "Mountain west",
                                                 "Appalachians")))))

Mypal <- c("#FF7F00","dodgerblue2","green4","gold1","#6A3D9A","gray70")


HUC2map = tm_shape(HUC2)+
  tm_fill(col = "Region", 
          palette = Mypal,
          alpha = 0.5)+
  tm_borders("grey10")+
  tm_shape(US_Bounds)+
  tm_borders("grey10")+
  tm_shape(HUC2)+
  tm_text("huc2", col = "black",
          size = 1.5, just = "right",bg.color = NA, remove.overlap = T)+
  tm_layout(frame = FALSE,
            legend.show = F)

# read in data
dat = fishdat

fish_tax <- read.csv("./data/FishTaxonomy_42024.csv")
fish_nat_stat <- read.csv("./data/NonnativeHUC2.csv",
                          colClasses = list("HUC2" = "character"))
gamey = read.csv("./Data/GameFishDesignation.csv")


# create an occurrence dataset for temporal analyses
dat_occur <- dat %>%
  dplyr::select(CollectionYear, HUC2, Agency,
                Luxilus.cornutus:Pteronotropis.metallicus) %>%
  group_by(CollectionYear, HUC2, Agency) %>%
  # calculate number of sites within regions and years
  mutate(n.sites = n()) %>%
  ungroup() %>%
  #sum presence/absence within region, years, and agency
  group_by(CollectionYear, HUC2, n.sites, Agency) %>%
  summarize_all(.funs = sum)%>%
  ungroup() %>%
  # create region-agency column
  mutate(reg_ag = paste(HUC2, Agency, sep = "_")) %>%
  dplyr::select(CollectionYear, HUC2, Agency, reg_ag, n.sites,
                Luxilus.cornutus:Pteronotropis.metallicus) 


##remove NAs for HUC2
dat_occur = dat_occur %>%
  filter(!is.na(HUC2)) %>%
  mutate(Agency = as.factor(Agency))


# make the dataset long in its orientation
dat_occur_long = dat_occur %>%
  pivot_longer(cols = Luxilus.cornutus:Pteronotropis.metallicus,
               names_to = "Species",
               values_to = "Occur") %>%
  #join families
  left_join(fish_tax %>% dplyr::select(Family, Species_period),
            join_by(Species == "Species_period")) %>%
  # make sure these columns are factors and not character
  mutate(across(c(HUC2, Agency,Species, Family), ~as.factor(.))) 
# filter out the few spp that don't have temp regimes associated

# 1) Evaluate families within HUC's to get at spatial differences

m_fam_gm <- bam(cbind(Occur, n.sites-Occur) ~ CollectionYear*Family*HUC2 +
                 Agency,
               data = dat_occur_long,
               family = "binomial")

summary(m_fam_gm)
anova(m_fam_gm)

##facetting in a specific order
##run this before your ggplot code
mygrid <- data.frame(
  code = c("01", "02", "03", "04", "05","06","07","08","09","10","11","12","13","14","15","16","17","18"),
  name = c("01", "02", "03", "04", "05","06","07","08","09","10","11","12","13","14","15","16","17","18"),
  row = c(1, 2, 3, 1, 2, 2, 1, 3, 1, 1, 2, 3, 3, 2, 3, 2, 1, 3),
  col = c(6, 6, 6, 5, 4, 5, 4, 5, 3, 2, 3, 4, 3, 2, 2, 1, 1, 1),
  stringsAsFactors = FALSE
)
mygrid = mygrid %>%
  mutate(code = paste("HUC2: ",code, sep = ""),
         name = paste("HUC2: ",name, sep = ""))

# extract data from the model for the plot
dat <- data.frame(emtrends(m_fam_gm, ~Family|HUC2,
                           var = "CollectionYear")) %>%
  mutate(CollectionYear.trend = CollectionYear.trend,
         lower.CL = lower.CL,
         upper.CL = upper.CL) %>%
  # filter out a few with huge errors
  filter(SE < 1) %>%
  # add column to work with geofacet
  mutate(name = HUC2,
         name = paste("HUC2: ",name, sep = "")) %>%
  group_by(HUC2)%>% 
  mutate(name2 = fct_reorder(Family, CollectionYear.trend)) %>%
  ungroup()%>%
  # drop families not significantly increasing/decreasing
  filter(lower.CL >0 | upper.CL <0) 


colorfacs = data.frame(name = unique(dat$HUC2)) %>%
  mutate(region = ifelse(name %in% c("02","03","08","12"),
                         "Coastal plains",
                         ifelse(name %in% c("04","05","07","09"),
                                "Midwest",
                                ifelse(name %in% c("13","14","15","16"),
                                       "Xeric",
                                       ifelse(name %in% c("10","11"),
                                              "Great plains",
                                              ifelse(name %in% c("17","18"),
                                                     "Mountain west",
                                                     "Appalachians"))))),
         CollectionYear.trend = 0,
         name2 = "name",
         name = paste("HUC2: ",name, sep = ""),
         HUC2 = name)



fig1 <- ggplot(data = dat, aes(y = name2, x = CollectionYear.trend))+
  geom_rect(data = colorfacs, inherit.aes = F,
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = region),alpha = 0.5) + 
  # facet_wrap(~HUC2, scales = "free_y")+
  facet_geo(~name, grid = mygrid, scales = "free_y")+
  
  #scale_x_continuous(trans = "log")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_errorbarh(aes(xmin = lower.CL,
                     xmax = upper.CL),
                 height = 0)+
  geom_point() +
  scale_fill_manual(values = c("#FF7F00","dodgerblue2","green4","gold1","#6A3D9A","gray70"),
                    name = "Region")+
  ylab("Family")+
  xlab("Trend (occurrence ~ time)")+
  theme_bw()+ 
  theme(axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size= 15),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15)) 

HUC2map_Grob = tmap_grob(HUC2map)

fulls5 = plot_grid(HUC2map_Grob,fig1,
                     nrow = 2,
                     labels = c("A","B"),
                     rel_heights = c(.85,1),
                     label_x = c(0.23,0))


ggsave("./figures/figFamHuc2.png", fulls5, dpi = 300, 
       height = 12, width = 15)


#temperature and thermal preference distirbutions####
streamTemps <- readRDS("./Data/powell-long-term-water-temperature-predictions.2024.08.08.rds")

# create dataset of mean stream temps for June and July
streamTemps <- streamTemps %>%
  dplyr::select(COMID, year, wt_pred) %>%
  group_by(COMID, year) %>%
  summarize(wt_pred = mean(wt_pred)) %>%
  ungroup() 

FD2 = streamTemps %>%
  mutate(YearG = ifelse(year < 1998 & year > 1992,
                        # year == 1990,
                        "Early",
                        ifelse(year > 2014,
                               # year == 2020,
                               "Late",
                               "Drop"
                        ))) %>%
  filter(YearG != "Drop") %>%
  mutate(PastRegime = ifelse(wt_pred < lowerbreak,
                             "Cold",
                             ifelse(wt_pred > upperbreak,
                                    "Warm",
                                    "Intermediate")))

FD2 %>% group_by(YearG) %>%
  summarize(m = mean(wt_pred),
            sd = sd(wt_pred))

##add means and SD to plots
predstreamtemps =  ggplot(FD2 ,aes(x = wt_pred))+#, fill = PastRegime))+
  annotate("rect", xmin = 0, xmax = lowerbreak, ymin = -.1, ymax = 1, 
           fill="blue", color=NA, alpha = 0.5)+
  annotate("rect", xmin = lowerbreak, xmax = upperbreak, ymin = -.1, ymax = 1, 
           fill="violet", color=NA, alpha = 0.5)+
  annotate("rect", xmin = upperbreak, xmax = 40, ymin = -.1, ymax = 1, 
           fill="red", color=NA, alpha = 0.5)+
  geom_density(trim = T,
               alpha = 0.25,
               adjust = 1,
               aes(linetype = YearG))+
  geom_errorbarh(data = FD2 %>% group_by(YearG) %>%
                   summarize(m = mean(wt_pred),
                             sd = sd(wt_pred)) %>%
                   ungroup() %>%
                   mutate(y = rep(c(.024,.012), each = 1)),
                 inherit.aes = F,
                 aes(y = y, group = YearG,
                     xmin = m - sd, xmax = m + sd),
                 height = 0)+
  geom_point(data = FD2 %>% group_by(YearG) %>%
               summarize(m = mean(wt_pred),
                         sd = sd(wt_pred)) %>%
               ungroup() %>%
               mutate(y = rep(c(.024,.012), each = 1)),
             inherit.aes = F,
             aes(x = m, y = y, shape = YearG),
             size = 3, fill = "white")+
  # scale_fill_manual(name = "Temperature regime",
  #                   values = c("blue","violet","red"),
  #                   guide= "none")+
  scale_x_continuous(labels = scales::label_number(suffix = "°C"))+
  scale_shape_manual(values = c(19,21),
                     labels = c("1993-1997", "2015-2019"),
                     name = "")+
  scale_linetype_manual(values = c("solid","dashed"),
                        labels = c("1993-1997", "2015-2019"),
                        name = "")+
  # scale_x_discrete(labels = c("1993-1997","2015-2019"))+
  xlab("Mean summer temperature")+
  ylab("Relative density")+
  coord_cartesian(ylim = c(.005,.15),
                  xlim = c(5,32))+
  theme_bw()+
  theme(legend.position = "none",
        legend.margin = margin(c(-5,0,0,0)),
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9),
        legend.title = element_text(size = 1),
        legend.background = element_rect(fill = "white", color = "black"))

###Abundance-based CWTP 

CWTP <- fishdatCPUE %>%
  dplyr::select(-wt_pred) %>%
  # join stream temps to fish data with temporally matched years
  left_join(streamTemps  %>%
              mutate(COMID = as.character(COMID)) %>%
              filter(COMID %in% fishdatCPUE$COMID),
            by = join_by(COMID == COMID, CollectionYear == year)) %>%
  mutate(PastRegime = ifelse(wt_pred < lowerbreak,
                             "Cold",
                             ifelse(wt_pred > upperbreak,
                                    "Warm",
                                    "Intermediate"))) %>%
  dplyr::select(SiteNumber, CollectionYear, HUC2, Agency, wt_pred, PastRegime,
                Luxilus.cornutus:Pteronotropis.metallicus) %>%
  pivot_longer(cols=Luxilus.cornutus:Pteronotropis.metallicus,
               values_to = "Obs",
               names_to = "species") %>%
  filter(Obs != 0) %>%
  left_join(FishLifeTraits %>%
              dplyr::select(Taxa, temperature),
            by = join_by(species == Taxa)) %>%
  mutate(ObsTherm = Obs * temperature) %>%
  group_by(SiteNumber, CollectionYear, HUC2, Agency, wt_pred, PastRegime) %>%
  summarize(CWTP  = weighted.mean(temperature, Obs),
            ObsS = sum(Obs),
            CWTP_a = CWTP / ObsS)


CWTP = CWTP %>%
  mutate(YearG = ifelse(CollectionYear < 1998,
                        "Early",
                        ifelse(CollectionYear > 2014,
                               "Late",
                               "Drop"
                        ))) %>%
  filter(YearG != "Drop")

CWTP %>% group_by(YearG) %>%
  summarize(m = mean(CWTP, na.rm = T),
            sd = sd(CWTP, na.rm = T))

##change the later mean point to open circle

##add means and SD to plots
CWTP_fig_a =   ggplot(CWTP,aes(x = CWTP))+#, fill = PastRegime))+
  geom_density(trim = T,
               alpha = 0.25,
               aes(linetype = YearG))+
  geom_errorbarh(data = CWTP %>% group_by(YearG) %>%
                   summarize(m = mean(CWTP, na.rm = T),
                             sd = sd(CWTP, na.rm = T)) %>%
                   ungroup() %>%
                   mutate(y = rep(c(.032,.016), each = 1)),
                 inherit.aes = F,
                 aes(y = y, group = YearG,
                     xmin = m - sd, xmax = m + sd),
                 height = 0)+
  geom_point(data = CWTP %>% group_by(YearG) %>%
               summarize(m = mean(CWTP, na.rm = T),
                         sd = sd(CWTP, na.rm = T)) %>%
               ungroup() %>%
               mutate(y = rep(c(.032,.016), each = 1)),
             inherit.aes = F,
             aes(x = m, y = y, shape = YearG),
             size = 3,
             fill = "white")+
  # scale_fill_manual(name = "Temperature regime",
  #                   values = c("blue","violet","red"),
  #                   guide= "none")+
  scale_x_continuous(labels = scales::label_number(suffix = "°C"))+
  scale_y_continuous(limits = c(-.005,0.2),
                     expand = c(0,0))+
  scale_shape_manual(values = c(19,21),
                     labels = c("1993-1997", "2015-2019"),
                     name = "")+
  scale_linetype_manual(values = c("solid","dashed"),
                        labels = c("1993-1997", "2015-2019"),
                        name = "")+
  xlab("Community-weighted thermal\npreference (abundance based)")+
  ylab("Relative density")+
  theme_bw()+
  theme(legend.position = "none",
        legend.margin = margin(c(-5,0,0,0)),
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9))

###Occurrence-based CWTP 

CWTP2 <- fishdat %>%
  # dplyr::select(-wt_pred) %>%
  # join stream temps to fish data with temporally matched years
  left_join(streamTemps  %>%
              mutate(COMID = as.character(COMID)) %>%
              filter(COMID %in% fishdatCPUE$COMID),
            by = join_by(COMID == COMID, CollectionYear == year)) %>%
  mutate(PastRegime = ifelse(wt_pred < lowerbreak,
                             "Cold",
                             ifelse(wt_pred > upperbreak,
                                    "Warm",
                                    "Intermediate"))) %>%
  dplyr::select(SiteNumber, CollectionYear, HUC2, Agency, wt_pred, PastRegime,
                Luxilus.cornutus:Pteronotropis.metallicus) %>%
  pivot_longer(cols=Luxilus.cornutus:Pteronotropis.metallicus,
               values_to = "Obs",
               names_to = "species") %>%
  filter(Obs != 0) %>%
  left_join(FishLifeTraits %>%
              dplyr::select(Taxa, temperature),
            by = join_by(species == Taxa)) %>%
  mutate(ObsTherm = Obs * temperature) %>%
  group_by(SiteNumber, CollectionYear, HUC2, Agency, wt_pred, PastRegime) %>%
  summarize(CWTP  = weighted.mean(temperature, Obs),
            ObsS = sum(Obs),
            CWTP_a = CWTP / ObsS)

CWTP2 = CWTP2 %>%
  mutate(YearG = ifelse(CollectionYear < 1998,
                        "Early",
                        ifelse(CollectionYear > 2014,
                               "Late",
                               "Drop"
                        ))) %>%
  filter(YearG != "Drop")

CWTP2 %>% group_by(YearG) %>%
  summarize(m = mean(CWTP, na.rm = T),
            sd = sd(CWTP, na.rm = T))

CWTP_fig_o =  ggplot(CWTP2,aes(x = CWTP))+#, fill = PastRegime))+
  geom_density(trim = T,
               alpha = 0.25,
               aes(linetype = YearG))+
  geom_errorbarh(data = CWTP %>% group_by(YearG) %>%
                   summarize(m = mean(CWTP, na.rm = T),
                             sd = sd(CWTP, na.rm = T)) %>%
                   ungroup() %>%
                   mutate(y = rep(c(.06,.03), each = 1)),
                 inherit.aes = F,
                 aes(y = y, group = YearG,
                     xmin = m - sd, xmax = m + sd),
                 height = 0)+
  geom_point(data = CWTP2 %>% group_by(YearG) %>%
               summarize(m = mean(CWTP, na.rm = T),
                         sd = sd(CWTP, na.rm = T)) %>%
               ungroup() %>%
               mutate(y = rep(c(.06,.03), each = 1)),
             inherit.aes = F,
             aes(x = m, y = y, shape = YearG),
             size = 3,
             fill = "white")+
  scale_x_continuous(labels = scales::label_number(suffix = "°C"))+
  scale_y_continuous(limits = c(-.005,0.375),
                     expand = c(0,0))+
  scale_shape_manual(values = c(19,21),
                     labels = c("1993-1997", "2015-2019"),
                     name = "")+
  scale_linetype_manual(values = c("solid","dashed"),
                        labels = c("1993-1997", "2015-2019"),
                        name = "")+
  xlab("Community-weighted thermal\npreference (occurrence based)")+
  ylab("Relative density")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.margin = margin(c(-5,0,0,0)),
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9))

cowplot::plot_grid(predstreamtemps,
                   CWTP_fig_a,
                   CWTP_fig_o,
                   allign = "hv",
                   nrow = 3,
                   ncol = 1,
                   labels = c("A","B","C"),
                   rel_heights = c(.8,.85,1))  

ggsave("./Figures/SuppFig_TempsDistribution.jpg", dpi = 300,
       height = 6,
       width = 2.5,
       units = "in")  


#Surface temperatures over last 70 years----

LandTempsNL = read.csv("./Data/ContUSSummerTemps.csv")

n1 = LandTempsNL %>%
  pivot_longer(cols = July:August,
               names_to = "Month",
               values_to = "Temp") %>%
  filter(Year %in% c(1950:2020)) %>%
  group_by(Year) %>%
  summarize(AvgSummerTemp = (mean(Temp) - 32) * 5/9)

mean((n1 %>% filter(Year %in% 1950:2000))$AvgSummerTemp)
sd((n1 %>% filter(Year %in% 1950:2000))$AvgSummerTemp)
mean((n1 %>% filter(Year %in% 1990:1994))$AvgSummerTemp)
sd((n1 %>% filter(Year %in% 1990:1994))$AvgSummerTemp)

##no difference between the 50 year average and 5 year average, 0.5F cooler
t.test((n1 %>% filter(Year %in% 1950:2000))$AvgSummerTemp,
       (n1 %>% filter(Year %in% 1990:1994))$AvgSummerTemp)

