library(finsyncR); 
library(StreamCatTools); 
library(tidyverse);
library(sf); library(nhdplusTools); library(vegan); library(FishLife)
library(adespatial); 

####Initial data generation and cleaning####
##read in occurance dataset
fishdat <- getFishData(dataType = "occur", taxonLevel = "Species", hybrids = F,
                       sharedTaxa = T)

##read in CPUE dataset
fishdatCPUE <- getFishData(dataType = "abun",
                           standardize = "CPUE",
                           taxonLevel = "Species",
                           hybrids = F,
                           sharedTaxa = T)

#calculate Species richness and TotalCPUE for each dataset, respectively
fishdat$SppRichness = rowSums(fishdat[,24:412])
fishdatCPUE$TotalCPUE <- rowSums(fishdatCPUE[,24:405])

fishdat$Year = fishdat$CollectionYear - min(fishdat$CollectionYear)  + 1

fishdatCPUE$Year = fishdatCPUE$CollectionYear - min(fishdatCPUE$CollectionYear)  + 1

fishdat$StreamOrder = as.character(fishdat$StreamOrder)
fishdatCPUE$StreamOrder = as.character(fishdatCPUE$StreamOrder)

saveRDS(fishdat, "./Data/Raw_fishdat.RDS")
saveRDS(fishdatCPUE, "./Data/Raw_fishdatCPUE.RDS")

##
fishdat = readRDS("./Data/Raw_fishdat.RDS")
fishdatCPUE = readRDS("./Data/Raw_fishdatCPUE.RDS")

##filter to just samples collected via electroshocking
nrow(fishdat) #6321
fishdat = fishdat %>% filter(SampleMethod == "Shocking")
6321 - nrow(fishdat) #644 diff (5677)

##remove samples not in contiguous USA
fishdat = fishdat %>% filter(Longitude_dd > -125)
5677 - nrow(fishdat) #27 diff (5650)

##randomly select 1 sample for each stream reach per year (if applicable)
set.seed(1)
fishdat <- fishdat %>%
  group_by(SiteNumber, CollectionYear) %>%
  slice(1)

5650 - nrow(fishdat) #455 diff (5195)

#CPUE
##filter to just samples collected via electroshocking
nrow(fishdatCPUE) #6044
fishdatCPUE = fishdatCPUE %>% filter(SampleMethod == "Shocking")
6044 - nrow(fishdatCPUE) #693 diff (5351)

##remove samples not in contiguous USA
fishdatCPUE = fishdatCPUE %>% filter(Longitude_dd > -125)
5351 - nrow(fishdatCPUE) #24 diff (5326)

##randomly select 1 sample for each stream reach per year (if applicable)
set.seed(1)
fishdatCPUE <- fishdatCPUE %>%
  group_by(SiteNumber, CollectionYear) %>%
  slice(1)

5326 - nrow(fishdatCPUE) #329 diff (4997)



####Conductivity data####
##USGS generated from TADA R package
usgsfish <- readRDS("./Data/USGSOccur_Conductivity.rds")
usgsfish2 <- readRDS("./Data/USGSCPUE_Conductivity.rds")
##EPA data taken directly from EPA NRSA website
EPA_wholecond <- readRDS("./Data/EPASites_Conductivity.rds")

fishdat <- fishdat %>%
  left_join(EPA_wholecond %>%
              ungroup() %>%
              group_by(UID) %>%
              slice(1) %>%
              ungroup() %>%
              mutate(UID = as.character(UID)) %>%
              dplyr::select(UID, CONDUCTIVITY2) %>%
              filter(UID %in% fishdat$SampleID),
            by = join_by("SampleID" == "UID")) %>%
  mutate(CollectionDate = as.Date(CollectionDate)) %>%
  left_join(usgsfish %>%
              group_by(SiteNumber, CollectionDate) %>%
              slice(1) %>%
              ungroup() %>%
              dplyr::select(Agency:COMID, UseCond)) %>%
  mutate(WholeConductivity = ifelse(is.na(CONDUCTIVITY2),
                                    UseCond, 
                                    CONDUCTIVITY2))



fishdatCPUE <- fishdatCPUE %>%
  left_join(EPA_wholecond %>%
              ungroup() %>%
              group_by(UID) %>%
              slice(1) %>%
              ungroup() %>%
              mutate(UID = as.character(UID)) %>%
              dplyr::select(UID, CONDUCTIVITY2) %>%
              filter(UID %in% fishdat$SampleID),
            by = join_by("SampleID" == "UID")) %>%
  mutate(CollectionDate = as.Date(CollectionDate)) %>%
  left_join(usgsfish2 %>%
              group_by(SiteNumber, CollectionDate) %>%
              slice(1) %>%
              ungroup() %>%
              dplyr::select(Agency:COMID, UseCond)) %>%
  mutate(WholeConductivity = ifelse(is.na(CONDUCTIVITY2),
                                    UseCond, 
                                    CONDUCTIVITY2))


#####HUC2 and broader HUC (12 or 14 digit; reachcode)####

flowdat = readRDS("./Data/HUC2dat.RDS")
fishdat <- fishdat %>%
  left_join(flowdat %>% group_by(comid) %>%
              slice(1) %>% ungroup() %>%
              rename("COMID" = "comid"))

##CPUE
fishdatCPUE = fishdatCPUE %>%
  filter(!is.na(COMID))

flowdatCPUE <- readRDS("./Data/HUC2datCPUE.RDS")

fishdatCPUE <- fishdatCPUE %>%
  # dplyr::select(-HUC2, -HUC2Clust) %>%
  left_join(flowdatCPUE %>% group_by(comid) %>%
              dplyr::slice(1) %>% ungroup() %>%
              rename("COMID" = "comid"))


####Land use data####
dat = bind_rows(fishdat %>%
                  dplyr::select(SiteNumber, CollectionYear),
                fishdatCPUE %>%
                  dplyr::select(SiteNumber, CollectionYear)) %>%
  group_by(SiteNumber, CollectionYear) %>%
  slice(1) %>%
  ungroup()

dat1 = dat %>%
  filter(SiteNumber %in% (finsyncR:::.allsitesCOMID %>% filter(!is.na(COMID)))$SiteNumber)

NLCDdat = finsyncR::getNLCDData(data = dat1,
                                scale = "Cat",
                                group = TRUE)

fishdat = fishdat %>%
  # dplyr::select(-PctCrop_Cat, -PctOpn_Cat, -PctUrb_Cat, -PctFstWtr_cat,
  #               -PctFst_Cat, -PctWater_Cat, -TopNLCD, -Landuse) %>%
  left_join(NLCDdat %>%
              filter(SiteNumber %in% fishdat$SiteNumber)) %>%
  rowwise() %>%
  mutate(PctFstWtr_cat = PctFst_Cat + PctWater_Cat) %>%
  ungroup() %>%
  mutate(TopNLCD = pmax(PctCrop_Cat, PctOpn_Cat, PctUrb_Cat, PctFstWtr_cat),
         Landuse = ifelse(TopNLCD == PctCrop_Cat,
                          "Agriculture",
                          ifelse(TopNLCD == PctOpn_Cat,
                                 "Grassland",
                                 ifelse(TopNLCD == PctUrb_Cat,
                                        "Urban",
                                        "Forest/Wetland"))))

fishdatCPUE = fishdatCPUE %>%
  # dplyr::select(-PctCrop_Cat, -PctOpn_Cat, -PctUrb_Cat, -PctFstWtr_cat,
  #               -PctFst_Cat, -PctWater_Cat, -TopNLCD, -Landuse) %>%
  left_join(NLCDdat %>%
              filter(SiteNumber %in% fishdatCPUE$SiteNumber)) %>%
  rowwise() %>%
  mutate(PctFstWtr_cat = PctFst_Cat + PctWater_Cat) %>%
  ungroup() %>%
  mutate(TopNLCD = pmax(PctCrop_Cat, PctOpn_Cat, PctUrb_Cat, PctFstWtr_cat),
         Landuse = ifelse(TopNLCD == PctCrop_Cat,
                          "Agriculture",
                          ifelse(TopNLCD == PctOpn_Cat,
                                 "Grassland",
                                 ifelse(TopNLCD == PctUrb_Cat,
                                        "Urban",
                                        "Forest/Wetland"))))


####gather %sand and predicted wetted width####

##attach COMIDs to site-numbers
data = bind_rows(fishdat%>%
                   dplyr::select(COMID),
                 fishdatCPUE %>%
                   dplyr::select(COMID)) %>%
  filter(!is.na(COMID)) %>%
  group_by(COMID) %>%
  slice(1) %>%
  ungroup()

comid = paste(data$COMID, collapse = ",")

#sand
streamcat_sand = sc_get_data(metric = "sand", comid = comid, aoi = "catchment")

streamcat_sand = streamcat_sand %>%
  dplyr::select(COMID, SANDCAT) %>%
  mutate(across(SANDCAT, ~./ 100))

fishdat = fishdat %>%
  left_join(streamcat_sand %>% filter(COMID %in% fishdat$COMID))

fishdatCPUE = fishdatCPUE %>%
  left_join(streamcat_sand %>% filter(COMID %in% fishdatCPUE$COMID))

#wetted width
streamcat_wettedwidth = sc_get_data(metric = "wettedwidth", comid = comid)

fishdat = fishdat %>%
  mutate(COMID = as.character(COMID)) %>%
  left_join(streamcat_wettedwidth %>%
              mutate(COMID = as.character(COMID)) %>%
              dplyr::select(COMID, WETTEDWIDTH) %>%
            filter(COMID %in% fishdat$COMID),
            by = join_by(COMID == COMID)) %>%
  mutate(PredictedWettedWidth_m = ifelse(is.na(WETTEDWIDTH),
                                         PredictedWettedWidth_m,
                                         WETTEDWIDTH)) %>%
  dplyr::select(-WETTEDWIDTH)

fishdatCPUE = fishdatCPUE %>%
  mutate(COMID = as.character(COMID)) %>%
  left_join(streamcat_wettedwidth %>%
              mutate(COMID = as.character(COMID)) %>%
              dplyr::select(COMID, WETTEDWIDTH) %>%
              filter(COMID %in% fishdatCPUE$COMID),
            by = join_by(COMID == COMID)) %>%
  mutate(PredictedWettedWidth_m = ifelse(is.na(WETTEDWIDTH),
                                         PredictedWettedWidth_m,
                                         WETTEDWIDTH)) %>%
  dplyr::select(-WETTEDWIDTH)


###Stream Temps####
streamTemps <- readRDS("./Data/powell-long-term-water-temperature-predictions.2024.08.08.rds")

COMID_temp_temporal <- function(df) {
  
  lm(wt_pred ~ year, data = df)
  
}

StreamTempTrends = streamTemps %>%
  mutate(year = year - 1989) %>%
  dplyr::select(COMID, year, wt_pred) %>%
  group_by(COMID, year) %>%
  summarize(wt_pred = mean(wt_pred)) %>%
  ungroup() %>%
  group_by(COMID) %>%
  nest() %>%
  mutate(model = map(data, COMID_temp_temporal)) %>%
  mutate(coefs = map(model,broom::tidy)) %>%
  unnest(coefs) %>%
  filter(term == "year") %>%
  dplyr::select(COMID,estimate, p.value)

SummerstreamTemps <- streamTemps %>%
  filter(year %in% c(1990:1994)) %>%
  dplyr::select(COMID, year, wt_pred) %>%
  group_by(COMID, year) %>%
  summarize(wt_pred = mean(wt_pred)) %>%
  ungroup() %>%
  group_by(COMID) %>%
  summarize(wt_pred = mean(wt_pred))

SummerstreamTemps = SummerstreamTemps %>% 
  left_join(StreamTempTrends)

# saveRDS(SummerstreamTemps, "./Data/SummerstreamTemps.RDS")

fishdat = fishdat %>%
  mutate(COMID = as.character(COMID)) %>%
  left_join(SummerstreamTemps %>%
              rename(wt_pred_new = wt_pred) %>%
              dplyr::select(COMID, wt_pred_new))

fishdatCPUE = fishdatCPUE %>%
  mutate(COMID = as.character(COMID)) %>%
  left_join(SummerstreamTemps %>%
              rename(wt_pred_new = wt_pred) %>%
              dplyr::select(COMID, wt_pred_new))


####Final data cleaning####
nrow(fishdat) #5195
fishdat  <- fishdat %>%
  mutate(PredictedWettedWidth_km = PredictedWettedWidth_m / 1000) %>%
  filter(complete.cases(WholeConductivity),
         complete.cases(PredictedWettedWidth_m),
         PredictedWettedWidth_m != 0,
         complete.cases(SANDCAT),
         complete.cases(Landuse),
         complete.cases(wt_pred_new),
         complete.cases(HUC2),
         StreamOrder != "0")

5195 - nrow(fishdat) #699 diff (4496)
  
fishdat  <- fishdat %>%
  filter(WholeConductivity < 20000)

4496 - nrow(fishdat) #1 diff (4495)

##remove temperature outliers
fishdat = fishdat %>% filter(wt_pred_new > 7)

4495 - nrow(fishdat) #4 diff (4491)

nrow(fishdatCPUE) #4996 (1 sample lost due to missing COMID)
fishdatCPUE <- fishdatCPUE  %>%
  mutate(PredictedWettedWidth_km = PredictedWettedWidth_m / 1000) %>%
  mutate(StreamOrder = as.character(StreamOrder)) %>%
  filter(complete.cases(WholeConductivity),
         complete.cases(PredictedWettedWidth_m),
         PredictedWettedWidth_m != 0,
         complete.cases(SANDCAT),
         complete.cases(Landuse),
         complete.cases(wt_pred_new),
         complete.cases(HUC2),
         StreamOrder != "0")
4996 - nrow(fishdatCPUE) #648 diff (4348)


fishdatCPUE <- fishdatCPUE  %>%
  filter(WholeConductivity < 20000)

4348 - nrow(fishdatCPUE) #1 diff (4347)

##filter based on sampling effort
fishdatCPUE = fishdatCPUE %>% filter(!is.infinite(TotalCPUE),
                       MethodEffort <= 150,
                       MethodEffort >= 3,
                       ReachLengthFished_m >= 10,
                       ReachLengthFished_m <=4000) 

4347 - nrow(fishdatCPUE) #104 diff (4243)

##remove temperature outliers
fishdatCPUE = fishdatCPUE %>% filter(wt_pred_new > 7)

4243 - nrow(fishdatCPUE) #3 diff (4240)

##make a categorical "year" variable to be used as categorical random effect
fishdat$YearC = as.character(fishdat$Year)
fishdatCPUE$YearC = as.character(fishdatCPUE$Year)

##convert TotalCPUE per m to TotalCPUE per 100m
fishdatCPUE$EffortOffset = (fishdatCPUE$MethodEffort * fishdatCPUE$ReachLengthFished_m)
fishdatCPUE$TotalCount = round(fishdatCPUE$TotalCPUE * (fishdatCPUE$EffortOffset),0)
fishdatCPUE$TotalCPUE100m = fishdatCPUE$TotalCPUE * 100
fishdatCPUE$ReachLengthFished_100m = fishdatCPUE$ReachLengthFished_m / 100

length(unique(fishdat$SiteNumber))
length(unique(fishdatCPUE$SiteNumber))

saveRDS(fishdat, "./Data/CleanedOccurrenceData.RDS")
saveRDS(fishdatCPUE, "./Data/CleanedAbundanceData.RDS")

####read in final filtered datasets####
fishdat <- readRDS("./Data/CleanedOccurrenceData.RDS")
fishdatCPUE <- readRDS("./Data/CleanedAbundanceData.RDS")
SummerstreamTemps <- readRDS("./Data/SummerstreamTemps.RDS")

fishdatCPUE = fishdatCPUE %>% left_join(SummerstreamTemps %>% dplyr::select(COMID, wt_pred,estimate))

####generate rarefied richness dataset####

fishdatCPUE$EffortOffset = (fishdatCPUE$MethodEffort * fishdatCPUE$ReachLengthFished_m)

fishdatCPUE$TotalCount = round(fishdatCPUE$TotalCPUE * fishdatCPUE$EffortOffset,0)

fishdatCPUE_spp = fishdatCPUE %>%
  filter(!is.na(TotalCount),
         is.finite(TotalCount),
         TotalCount > 24) %>%
  dplyr::select(Luxilus.cornutus:Pteronotropis.metallicus, EffortOffset) %>%
  mutate(across(Luxilus.cornutus:Pteronotropis.metallicus, ~round(. * EffortOffset, 0))) %>%
  dplyr::select(-EffortOffset)

S <- specnumber(fishdatCPUE_spp) # observed number of species
(raremax <- min(rowSums(fishdatCPUE_spp)))
set.seed(1)
Srare <- rarefy(fishdatCPUE_spp, raremax)

fishdat_rarefiedRich = (fishdatCPUE %>%
                          filter(!is.na(TotalCount),
                                 is.finite(TotalCount),
                                 TotalCount > 24))

fishdat_rarefiedRich$RarefiedRichness = Srare

####Generate LCBD values ####

fishdat$Richness = rowSums(fishdat[,24:412])

fishdat_LCBD = fishdat %>%
  filter(Richness > 0)

h1 = h2 = fishcomhld = data.frame()

uHUC2 = unique(fishdat$HUC2)

for(i in 1:18){
  h1 <- fishdat_LCBD %>%
    filter(HUC2 == uHUC2[i])
  fishcomhld = data.frame(h1[,24:412])
  fishcomhld = fishcomhld[,colSums(fishcomhld) > 0]
  fishdist1 <- vegdist(data.frame(h1[,24:412]), dist = "raup")
  
  h1$LCBD = LCBD.comp(fishdist1, sqrt.D = T)$LCBD
  ##center and scale (z-score) LCBD within regions
  h1$LCBD_s = as.numeric(scale(h1$LCBD) )
  h2 = rbind(h2,h1)
}

fishdat_LCBD = fishdat_LCBD %>%
  left_join(h2 %>% dplyr::select(SiteNumber, CollectionDate, LCBD,LCBD_s))

##rescale to between 0 and 1, for ease of interpretation
fishdat_LCBD$LCBD_s2 = scales::rescale(fishdat_LCBD$LCBD_s, to = c(0,1))

####Gather functional diversity data####

SESBoth = readRDS("./Data/SESWholeCommFinal.RDS")

#functional dispersion and Standard effect size data
fishdat <- fishdat %>%
  left_join(SESBoth %>% filter(SiteNumber %in% fishdat$SiteNumber)  %>%
              group_by(SiteNumber, CollectionYear) %>%
              slice(1) %>%
              dplyr::select(SiteNumber, CollectionYear, SES, fdis),
            by = join_by(SiteNumber == SiteNumber,
                         CollectionYear == CollectionYear))


###native, non-game####
gamey = read.csv("./Data/GameFishDesignation.csv")
nonnative1 = read.csv("./Data/USGSinput2_update_screened_nonnative.csv")
nonnative2 = read.csv("./Data/USGSinput1_update_screened_nonnative.csv")
nonnative = rbind(nonnative1, nonnative2)

nonnative = nonnative %>%
  mutate(HUC8 = as.character(str_pad(HUC8, 8, pad = "0")),
         Scientific.Name = gsub(" ","\\.",Scientific.Name),
         HUC8_Species = paste(HUC8, Scientific.Name, sep = "_")) %>%
  filter(Native == "No")

gamey$Species = gsub(" ", ".", gamey$Species)


fishdat_native = fishdat %>%
  mutate(reachcode = as.character(str_pad(reachcode, 14, pad = "0")),
         HUC8 = substr(reachcode,1,8)) %>%
  pivot_longer(cols = Luxilus.cornutus:Pteronotropis.metallicus,
               names_to = "Species",
               values_to = "Occur") %>%
  mutate(HUC8_Species = paste(HUC8, Species, sep = "_")) %>%
  filter(!(HUC8_Species %in% nonnative$HUC8_Species)) %>%
  filter(!(Species %in% gamey$Species)) %>%
  filter(Species != "Cyprinus.carpio") %>%
  dplyr::select(-HUC8_Species) %>%
  pivot_wider(names_from = Species,
              values_from = Occur,
              values_fill = 0)


fishdatCPUE_native = fishdatCPUE %>%
  mutate(reachcode = as.character(str_pad(reachcode, 14, pad = "0")),
         HUC8 = substr(reachcode,1,8)) %>%
  pivot_longer(cols = Luxilus.cornutus:Pteronotropis.metallicus,
               names_to = "Species",
               values_to = "Occur") %>%
  mutate(HUC8_Species = paste(HUC8, Species, sep = "_")) %>%
  filter(!(HUC8_Species %in% nonnative$HUC8_Species)) %>%
  filter(!(Species %in% gamey$Species)) %>%
  filter(Species != "Cyprinus.carpio") %>%
  dplyr::select(-HUC8_Species) %>%
  pivot_wider(names_from = Species,
              values_from = Occur,
              values_fill = 0)

fishdatCPUE_native$TotalCPUE = rowSums(fishdatCPUE_native[,50:379])
fishdatCPUE_native$TotalCPUE100m = fishdatCPUE_native$TotalCPUE * 100
fishdatCPUE_native$ReachLengthFished_100m = fishdatCPUE_native$ReachLengthFished_m / 100
fishdatCPUE_native$YearC = as.character(fishdatCPUE_native$Year)



fishdatCPUE_native$EffortOffset = (fishdatCPUE_native$MethodEffort * fishdatCPUE_native$ReachLengthFished_m)

fishdatCPUE_native$TotalCount = round(fishdatCPUE_native$TotalCPUE * fishdatCPUE_native$EffortOffset,0)

fishdatCPUE_spp_nng = fishdatCPUE_native %>%
  filter(!is.na(TotalCount),
         is.finite(TotalCount),
         TotalCount > 24) %>%
  dplyr::select(Luxilus.cornutus:Pteronotropis.metallicus, EffortOffset) %>%
  mutate(across(Luxilus.cornutus:Pteronotropis.metallicus, ~round(. * EffortOffset, 0))) %>%
  dplyr::select(-EffortOffset)

S <- specnumber(fishdatCPUE_spp_nng) # observed number of species
(raremax <- min(rowSums(fishdatCPUE_spp_nng)))
set.seed(1)
Srare <- rarefy(fishdatCPUE_spp_nng, raremax)

fishdat_rarefiedRich_nng = (fishdatCPUE_native %>%
                              filter(!is.na(TotalCount),
                                     is.finite(TotalCount),
                                     TotalCount > 24))

fishdat_rarefiedRich_nng$RarefiedRichness = Srare


#####functional diversity ####
SESBoth_nng = readRDS("./Data/SESNNGFinal.RDS")


fishdat_native <- fishdat_native %>%
  dplyr::select(-SES, -fdis) %>%
  left_join(SESBoth_nng %>% filter(SiteNumber %in% fishdat_native$SiteNumber)  %>%
              group_by(SiteNumber, CollectionYear) %>%
              slice(1) %>%
              dplyr::select(SiteNumber, CollectionYear, SES, fdis),
            by = join_by(SiteNumber == SiteNumber,
                         CollectionYear == CollectionYear))


#####beta diversity####

fishdat_native$Richness = rowSums(fishdat_native[,45:380])

fishdat_native_LCBD = fishdat_native %>%
  filter(Richness > 0)

h1 = h2 = data.frame()

uHUC2 = unique(fishdat_native_LCBD$HUC2)

for(i in 1:18){
  h1 <- fishdat_native_LCBD %>%
    filter(HUC2 == uHUC2[i])
  
  fishdist1 <- vegdist(data.frame(h1[,45:380]), dist = "raup")
  
  h1$LCBD2 = LCBD.comp(fishdist1, sqrt.D = T)$LCBD
  h1$LCBD2_s = as.numeric(scale(h1$LCBD2) )
  
  h2 = rbind(h2,h1)
  
}

fishdat_native_LCBD = fishdat_native_LCBD %>%
  dplyr::select(-any_of(c("LCBD2", "LCBD2_s"))) %>%
  left_join(h2 %>% dplyr::select(SiteNumber, CollectionDate, LCBD2,LCBD2_s))

fishdat_native_LCBD$LCBD2_s2 = scales::rescale(fishdat_native_LCBD$LCBD2_s, to = c(0,1))


###non-native, game####
fishdat_notnative = fishdat %>%
  mutate(reachcode = as.character(str_pad(reachcode, 14, pad = "0")),
         HUC8 = substr(reachcode,1,8)) %>%
  pivot_longer(cols = Luxilus.cornutus:Pteronotropis.metallicus,
               names_to = "Species",
               values_to = "Occur") %>%
  mutate(HUC8_Species = paste(HUC8, Species, sep = "_")) %>%
  filter((HUC8_Species %in% nonnative$HUC8_Species) |
           (Species %in% gamey$Species) |
           Species == "Cyprinus.carpio") %>%
  dplyr::select(-HUC8_Species) %>%
  pivot_wider(names_from = Species,
              values_from = Occur,
              values_fill = 0)


fishdatCPUE_notnative = fishdatCPUE %>%
  mutate(reachcode = as.character(str_pad(reachcode, 14, pad = "0")),
         HUC8 = substr(reachcode,1,8)) %>%
  pivot_longer(cols = Luxilus.cornutus:Pteronotropis.metallicus,
               names_to = "Species",
               values_to = "Occur") %>%
  mutate(HUC8_Species = paste(HUC8, Species, sep = "_")) %>%
  filter((HUC8_Species %in% nonnative$HUC8_Species) |
           (Species %in% gamey$Species) |
           Species == "Cyprinus.carpio") %>%
  dplyr::select(-HUC8_Species) %>%
  pivot_wider(names_from = Species,
              values_from = Occur,
              values_fill = 0)


fishdatCPUE_notnative$TotalCPUE = rowSums(fishdatCPUE_notnative[,50:347])
fishdatCPUE_notnative$TotalCPUE100m = fishdatCPUE_notnative$TotalCPUE * 100
fishdatCPUE_notnative$ReachLengthFished_100m = fishdatCPUE_notnative$ReachLengthFished_m / 100
fishdatCPUE_notnative$EffortOffset = (fishdatCPUE_notnative$MethodEffort * fishdatCPUE_notnative$ReachLengthFished_m)
fishdatCPUE_notnative$TotalCount = round(fishdatCPUE_notnative$TotalCPUE * fishdatCPUE_notnative$EffortOffset,0)

fishdatCPUE_spp_notnative = fishdatCPUE_notnative %>%
  filter(!is.na(TotalCount),
         is.finite(TotalCount),
         TotalCount > 24) %>%
  dplyr::select(Micropterus.dolomieu:Cottus.beldingii, EffortOffset) %>%
  mutate(across(Micropterus.dolomieu:Cottus.beldingii, ~round(. * EffortOffset, 0))) %>%
  dplyr::select(-EffortOffset)

S <- specnumber(fishdatCPUE_spp_notnative) # observed number of species
(raremax <- min(rowSums(fishdatCPUE_spp_notnative)))
set.seed(1)
Srare <- rarefy(fishdatCPUE_spp_notnative, raremax)

fishdat_rarefiedRich_notnative = (fishdatCPUE_notnative %>%
                                    filter(!is.na(TotalCount),
                                           is.finite(TotalCount),
                                           TotalCount > 24))

fishdat_rarefiedRich_notnative$RarefiedRichness = Srare

#####functional diversity####
SESBoth_notnative = readRDS("./Data/SESgame_notnativeFinal.RDS")

fishdat_notnativeSES <- fishdat_notnative %>%
  dplyr::select(-SES, -fdis) %>%
  left_join(SESBoth_notnative %>% filter(SiteNumber %in% fishdat_notnative$SiteNumber)  %>%
              group_by(SiteNumber, CollectionYear) %>%
              slice(1) %>%
              dplyr::select(SiteNumber, CollectionYear, SES, fdis),
            by = join_by(SiteNumber == SiteNumber,
                         CollectionYear == CollectionYear))

#####beta diversity####

fishdat_notnative$Richness = rowSums(fishdat_notnative[,47:350])

fishdat_notnative_LCBD = fishdat_notnative %>%
  filter(Richness > 0)


h1 = h2 = data.frame()

uHUC2 = unique(fishdat_notnative_LCBD$HUC2)

for(i in 1:18){
  h1 <- fishdat_notnative_LCBD %>%
    filter(HUC2 == uHUC2[i])
  
  fishdist1 <- vegdist(data.frame(h1[,47:350]), dist = "raup")
  
  h1$LCBD2 = LCBD.comp(fishdist1, sqrt.D = T)$LCBD
  h1$LCBD2_s = as.numeric(scale(h1$LCBD2) )
  
  h2 = rbind(h2,h1)
  
}

fishdat_notnative_LCBD = fishdat_notnative_LCBD %>%
  dplyr::select(-any_of(c("LCBD2", "LCBD2_s"))) %>%
  left_join(h2 %>% dplyr::select(SiteNumber, CollectionDate, LCBD2,LCBD2_s))

fishdat_notnative_LCBD$LCBD2_s2 = scales::rescale(fishdat_notnative_LCBD$LCBD2_s, to = c(0,1))

###Fish traits####

data(FishBase_and_Morphometrics)


edge_names = c( FishBase_and_Morphometrics$tree$tip.label,
                FishBase_and_Morphometrics$tree$node.label[-1] ) # Removing root

#
specieslist = c(colnames(fishdat)[24:412],"Etheostoma tetrazonum",
                'Hypophthalmichthys nobilis')
specieslist = gsub("\\.", " ",specieslist)

##species names updates from Fishbase
#Lampetra appendix -> Lethenteron appendix
#Notropis amplamala -> Ericymba amplamala
#Etheostoma chlorosoma ->	Etheostoma chlorosomum
##several species not found in the dataset, so getting Genus-level traits for those
specieslist[which(specieslist == "Lampetra appendix")] = "Lethenteron appendix"
specieslist[which(specieslist == "Notropis amplamala")] = "Ericymba amplamala"
specieslist[which(specieslist == "Etheostoma chlorosoma")] = "Etheostoma chlorosomum"
specieslist = c(specieslist,"Moxostoma", "Notropis",  "Hybopsis",  "Cottus")

which_g = 0

for(i in 1:length(specieslist)){
  which_g[i] = match( specieslist[i], edge_names )
}

fulldat <- data.frame(FishBase_and_Morphometrics$beta_gv)
fulldat$Taxa = row.names(fulldat)
fulldat$Taxa  = gsub("\\.", " ",fulldat$Taxa )

FishLifeTraits = fulldat %>%
  filter(Taxa %in% specieslist) %>%
  mutate(freq = ifelse(Taxa %in% c("Hybopsis","Moxostoma"),
                       2,
                       ifelse(Taxa == "Notropis",
                              3,
                              1))) %>%
  uncount(freq)

##need to be sure to join the 8? species to their genera. Should be easy enough.
##Should be able to just to an "anti join" to filter to those taxa in the species list
##that are not in the FishLifeTraits list, then just match the species names to the genus names
#
specieslist[is.na(which_g)]

##these are the species no found in the dataset
FishLifeTraits$Taxa[384:391] = c("Cottus chattahoochee",
                             "Moxostoma duquesnei",
                             "Moxostoma collapsum",
                             "Hybopsis winchelli",
                             "Hybopsis rubrifrons",
                             "Notropis percobromus",
                             "Notropis candidus",
                             "Notropis alborus")
rownames(FishLifeTraits) = gsub(" ",".",FishLifeTraits$Taxa)

##change "base" to actual values
colnames(FishLifeTraits)[c(18,21,26,29)] = c("spawning_typenonguarders",
                                         "habitatdemersal",
                                         "feeding_modegeneralist",
                                         "body_fusiform_normal")

colnames(FishLifeTraits)

FishLifeTraits = FishLifeTraits %>%
  dplyr::select(log.age_max.,
                trophic_level,
                temperature, 
                any_of(contains("spawning")),
                log.length_max.,
                log.age_maturity.,
                log.growth_coefficient.,
                log.fecundity.)

rownames(FishLifeTraits)[match(c("Lethenteron.appendix",
                             "Ericymba.amplamala",
                             "Etheostoma.chlorosomum"),
                           rownames(FishLifeTraits))] =
  c("Lampetra.appendix",
    "Notropis.amplamala",
    "Etheostoma.chlorosoma")
FishLifeTraits$Taxa = row.names(FishLifeTraits)


beta_iv = FishBase_and_Morphometrics$beta_gv[ FishBase_and_Morphometrics$g_i, ]

set.seed(6132014)
a <- archetypes::archetypes(beta_iv, 3, verbose = TRUE)

OPE = data.frame(a$alphas, species = row.names(beta_iv)) %>%
  group_by(species) %>%
  slice(1)

colnames(OPE)[1:3] = c("Equilibrium","Periodic","Opportunist")
# View(OPE[grepl("Scaphirhynchus",OPE$species ),])


##missing a handful of species; need to make some linkages
OPE_sub = OPE %>%
  filter(species %in% c(specieslist))


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


OPE_subfin %>%
  filter(grepl("piceus", species))
  filter(species %in% c("Ctenopharyngodon idella",
                        "Mylopharyngodon piceus"))

###
#prop non-native/game

fishdatCPUE_native = fishdatCPUE_native %>%
  left_join(
    fishdatCPUE_notnative %>%
      dplyr::select(SiteNumber, CollectionDate, TotalCPUE) %>%
      rename(NNaG_TotCPUE = TotalCPUE) %>%
      left_join(fishdatCPUE %>%
                  dplyr::select(SiteNumber, CollectionDate, TotalCPUE)) %>%
      mutate(PropNNaG = (NNaG_TotCPUE / TotalCPUE)) %>%
      dplyr::select(SiteNumber, CollectionDate, PropNNaG)
    )


fishdat_native = fishdat_native %>%
  left_join(
    fishdatCPUE_notnative %>%
      dplyr::select(SiteNumber, CollectionDate, TotalCPUE) %>%
      rename(NNaG_TotCPUE = TotalCPUE) %>%
      left_join(fishdatCPUE %>%
                  dplyr::select(SiteNumber, CollectionDate, TotalCPUE)) %>%
      mutate(PropNNaG = (NNaG_TotCPUE / TotalCPUE)) %>%
      dplyr::select(SiteNumber, CollectionDate, PropNNaG)
  )

fishdat_native_LCBD = fishdat_native_LCBD %>%
  left_join(
    fishdatCPUE_notnative %>%
      dplyr::select(SiteNumber, CollectionDate, TotalCPUE) %>%
      rename(NNaG_TotCPUE = TotalCPUE) %>%
      left_join(fishdatCPUE %>%
                  dplyr::select(SiteNumber, CollectionDate, TotalCPUE)) %>%
      mutate(PropNNaG = (NNaG_TotCPUE / TotalCPUE)) %>%
      dplyr::select(SiteNumber, CollectionDate, PropNNaG)
  )

fishdat_rarefiedRich_nng = fishdat_rarefiedRich_nng %>%
  left_join(
    fishdatCPUE_notnative %>%
      dplyr::select(SiteNumber, CollectionDate, TotalCPUE) %>%
      rename(NNaG_TotCPUE = TotalCPUE) %>%
      left_join(fishdatCPUE %>%
                  dplyr::select(SiteNumber, CollectionDate, TotalCPUE)) %>%
      mutate(PropNNaG = (NNaG_TotCPUE / TotalCPUE)) %>%
      dplyr::select(SiteNumber, CollectionDate, PropNNaG)
  )

##whole community occurrence data, probably

#OPE for each
##pivot wholecomm longer, then left_join to OPE_sub, filter out occur == 0,
##then group_by everything and take sum of each, get n(), divide sums by n()

fishdatOPE = fishdat %>%
  pivot_longer(cols = Luxilus.cornutus:Pteronotropis.metallicus,
               names_to = "Species",
               values_to = "Occur") %>%
  left_join(OPE_subfin %>% mutate(species = gsub(" ", ".", species)),
            by = join_by(Species == species)) %>%
  filter(Occur == 1) %>%
  group_by(across(c(-Species,-Occur, -Equilibrium, -Periodic, -Opportunist))) %>%
  summarize(Equilibrium = sum(Equilibrium, na.rm = T),
            Periodic = sum(Periodic, na.rm = T),
            Opportunist = sum(Opportunist, na.rm = T),
            nspecies = n()) %>%
  mutate(Equilibrium = Equilibrium / nspecies,
         Periodic = Periodic / nspecies,
         Opportunist = Opportunist / nspecies)  %>%
  mutate(across(Equilibrium:Opportunist, ~((. * (389 - 1)) + 0.5) / 389))
  # pivot_longer(cols = Equilibrium:Opportunist,
  #              names_to = "LifeHistory",
  #              values_to = "Likelihood")


#Abundance based
fishdatOPE = fishdatCPUE %>%
  pivot_longer(cols = Luxilus.cornutus:Pteronotropis.metallicus,
               names_to = "Species",
               values_to = "Occur") %>%
  left_join(OPE_subfin %>% mutate(species = gsub(" ", ".", species)),
            by = join_by(Species == species)) %>%
  filter(Occur > 0) %>%
  group_by(across(c(-Species,-Occur, -Equilibrium, -Periodic, -Opportunist))) %>%
  summarize(Equilibrium = weighted.mean(Equilibrium, Occur),
            Periodic = weighted.mean(Periodic, Occur),
            Opportunist = weighted.mean(Opportunist, Occur)) %>%
  # mutate(Equilibrium = Equilibrium ,
  #        Periodic = Periodic / nspecies,
  #        Opportunist = Opportunist / nspecies)  %>%
  mutate(across(Equilibrium:Opportunist, ~((. * (4226 - 1)) + 0.5) / 4226))
  # pivot_longer(cols = Equilibrium:Opportunist,
  #              names_to = "LifeHistory",
  #              values_to = "Likelihood")


#save datasets----



HUC2pres = cbind(fishdat %>%
  group_by(HUC2) %>%
  summarize(across(Luxilus.cornutus:Pteronotropis.metallicus, ~sum(.))),
fishdatCPUE %>%
  group_by(HUC2) %>%
  summarize(Etheostoma.tetrazonum = sum(Etheostoma.tetrazonum),
            Hypophthalmichthys.nobilis = sum(Hypophthalmichthys.nobilis)) %>%
  ungroup() %>%
  dplyr::select(-HUC2)) %>%
  mutate(across(Luxilus.cornutus:Hypophthalmichthys.nobilis, ~ifelse(. > 0,
                                                                     1, 0))) %>%
  pivot_longer(cols = Luxilus.cornutus:Hypophthalmichthys.nobilis,
               values_to = "Pres",
               names_to = "Species")

HUC2presNN = HUC2pres %>% left_join(nonnative %>%
  mutate(HUC2 = substr(HUC8, 1,2)) %>%
  group_by(Scientific.Name, HUC2) %>%
  summarize(Native = paste(Native, collapse = '_')) %>%
  mutate(Native = 'Non-native') %>%
  ungroup() %>%
  rename(Species = Scientific.Name)) %>%
  mutate(Native = ifelse(is.na(Native) & Pres == 1,
                         "Native",
                         Native)) %>%
  dplyr::select(-Pres) %>%
  filter(!is.na(Native)) %>%
  pivot_wider(names_from = HUC2,
              values_from = Native,
              names_prefix = "HUC",
              values_fill = "")




fishdesignations = OPE_subfin %>%
  mutate(species = gsub(" ", ".", species)) %>%
  left_join(gamey, by = join_by("species" == "Species")) %>%
  rename(GameFish = Status,
         Species = species) %>%
  left_join(HUC2presNN) %>%
  mutate(GameFish = ifelse(is.na(GameFish),
                           "Non-game",
                           "Game"),
         Species = gsub("\\."," ", Species)) %>%
  dplyr::select(Species, Equilibrium:Opportunist, GameFish, HUC01:HUC18)


write.csv(fishdesignations, "./Data/fishdesignations.csv", row.names = F)



###
whol_passes = read.csv("./Data/whol_passes.csv")

fishdatCPUE = fishdatCPUE %>%
left_join(whol_passes)


fishdatCPUE = fishdatCPUE %>%
  mutate(Npass = ifelse(is.na(Npass)|Npass == 0,
                        1,
                        Npass))


fishdat_rarefiedRich = fishdat_rarefiedRich %>%
  left_join(whol_passes)%>%
  mutate(Npass = ifelse(is.na(Npass)|Npass == 0,
                        1,
                        Npass))

fishdat = fishdat %>%
  left_join(whol_passes)%>%
  mutate(Npass = ifelse(is.na(Npass)|Npass == 0,
                        1,
                        Npass))

fishdat_LCBD = fishdat_LCBD %>%
  left_join(whol_passes)%>%
  mutate(Npass = ifelse(is.na(Npass)|Npass == 0,
                        1,
                        Npass))

fishdatOPE = fishdatOPE %>%
  left_join(whol_passes)%>%
  mutate(Npass = ifelse(is.na(Npass)|Npass == 0,
                        1,
                        Npass))

fishdatCPUE_native = fishdatCPUE_native %>%
  left_join(whol_passes)%>%
  mutate(Npass = ifelse(is.na(Npass)|Npass == 0,
                        1,
                        Npass))

fishdat_rarefiedRich_nng = fishdat_rarefiedRich_nng %>%
  left_join(whol_passes)%>%
  mutate(Npass = ifelse(is.na(Npass)|Npass == 0,
                        1,
                        Npass))

fishdat_native = fishdat_native %>%
  left_join(whol_passes)%>%
  mutate(Npass = ifelse(is.na(Npass)|Npass == 0,
                        1,
                        Npass))

fishdat_native_LCBD = fishdat_native_LCBD %>%
  left_join(whol_passes)%>%
  mutate(Npass = ifelse(is.na(Npass)|Npass == 0,
                        1,
                        Npass))

fishdatCPUE_notnative = fishdatCPUE_notnative %>%
  left_join(whol_passes)%>%
  mutate(Npass = ifelse(is.na(Npass)|Npass == 0,
                        1,
                        Npass))

fishdat_rarefiedRich_notnative = fishdat_rarefiedRich_notnative %>%
  left_join(whol_passes)%>%
  mutate(Npass = ifelse(is.na(Npass)|Npass == 0,
                        1,
                        Npass))

fishdat_notnativeSES = fishdat_notnativeSES %>%
  left_join(whol_passes)%>%
  mutate(Npass = ifelse(is.na(Npass)|Npass == 0,
                        1,
                        Npass))

fishdat_notnative_LCBD = fishdat_notnative_LCBD %>%
  left_join(whol_passes)%>%
  mutate(Npass = ifelse(is.na(Npass)|Npass == 0,
                        1,
                        Npass))