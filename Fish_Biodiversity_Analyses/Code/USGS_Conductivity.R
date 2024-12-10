remotes::install_github("USEPA/TADA",
                        ref = "develop",
                        dependencies = TRUE
)

library(TADA)
library(tidyverse)


testUSGSsites <- unique((fishdat1 %>% filter(Agency == "USGS"))$SiteNumber)


testUSGSsites2 <- unique((fishdatCPUE %>% filter(Agency == "USGS"))$SiteNumber)

allUSGSsites = unique(c(testUSGSsites,testUSGSsites2))

Allconductance <- TADA_DataRetrieval(characteristicName = c("Specific conductance",
                                                            "Specific conductivity",
                                                            "Conductivity",
                                                            "Conductance"),
                                        startDate = "1993-01-01",
                                        endDate = "2020-01-01",
                                        sampleMedia = c("Water", "water"),
                                        siteType = "Stream",
                                     siteid = allUSGSsites,
                                      applyautoclean = F)


##need to do a temporal match of the data
##So, Allconductance ActivityEndDate or ActivityStartDate need to fall within 2 weeks of the observed data

library(lubridate)

usgsfish <- fishdat1 %>%
  filter(Agency == "USGS") %>%
  ungroup() %>%
  mutate(CollectionDate = as.Date(CollectionDate),
         # MonthEarlier = CollectionDate -14,
         # MonthLater = CollectionDate +14
         MonthEarlier = CollectionDate %m-% months(2),
         MonthLater = CollectionDate %m+% months(2)
  )

usgsfish2 <- fishdatCPUE %>%
  filter(Agency == "USGS") %>%
  ungroup() %>%
  mutate(CollectionDate = as.Date(CollectionDate),
         # MonthEarlier = CollectionDate -14,
         # MonthLater = CollectionDate +14
         MonthEarlier = CollectionDate %m-% months(2),
         MonthLater = CollectionDate %m+% months(2)
  )

allusgsfish <- bind_rows(usgsfish %>% dplyr::select(SiteNumber, CollectionDate, MonthEarlier, MonthLater),
          usgsfish2 %>% dplyr::select(SiteNumber, CollectionDate, MonthEarlier, MonthLater)
) %>%
  group_by(SiteNumber, CollectionDate, MonthEarlier, MonthLater) %>%
  slice(1)



Allconductance2 = Allconductance %>%
  mutate(ActivityEndDate = as.Date(ActivityEndDate),
         ActivityStartDate = as.Date(ActivityStartDate)) %>%
  left_join(allusgsfish %>% dplyr::select(SiteNumber, CollectionDate, MonthEarlier, MonthLater),
            by = join_by("MonitoringLocationIdentifier" == "SiteNumber"))

Allcondsub <- Allconductance2 %>%
  filter(ActivityStartDate > MonthEarlier & ActivityStartDate <MonthLater |
           ActivityEndDate > MonthEarlier & ActivityEndDate <MonthLater) %>%
  dplyr::select(-MonthEarlier, -MonthLater)

Allcondsub <- TADA_AutoClean(Allcondsub)

AllcondsubCleaned <- Allcondsub %>%
  dplyr::select(ActivityStartDate,MonitoringLocationIdentifier,
                CollectionDate, TADA.ComparableDataIdentifier, 
                TADA.ResultMeasureValue, TADA.ResultMeasureValueDataTypes.Flag,
                ResultValueTypeName)

MeanConductivity <- AllcondsubCleaned %>%
  mutate(daydiff = as.numeric(ActivityStartDate - CollectionDate)) %>%
  group_by(MonitoringLocationIdentifier,CollectionDate) %>%
  mutate(minday = min(daydiff)) %>%
  ungroup() %>%
  filter(daydiff == minday) %>%
  group_by(MonitoringLocationIdentifier,CollectionDate) %>%
  summarize(Cond = mean(TADA.ResultMeasureValue, na.rm = T))


usgsfish <- usgsfish %>%
  mutate(SiteNumberCollectionDate = paste(SiteNumber, CollectionDate, sep = "_"))

usgsfish <- usgsfish %>%
  left_join(MeanConductivity %>% 
              ungroup() %>%
              mutate(SiteNumberCollectionDate = paste(MonitoringLocationIdentifier, CollectionDate, sep = "_")) %>%
              dplyr::select(-MonitoringLocationIdentifier, -CollectionDate) %>%
              filter(SiteNumberCollectionDate %in% usgsfish$SiteNumberCollectionDate)
              )


usgsfish2 <- usgsfish2 %>%
  mutate(SiteNumberCollectionDate = paste(SiteNumber, CollectionDate, sep = "_"))

usgsfish2 <- usgsfish2 %>%
  left_join(MeanConductivity %>% 
              ungroup() %>%
              mutate(SiteNumberCollectionDate = paste(MonitoringLocationIdentifier, CollectionDate, sep = "_")) %>%
              dplyr::select(-MonitoringLocationIdentifier, -CollectionDate) %>%
              filter(SiteNumberCollectionDate %in% usgsfish2$SiteNumberCollectionDate)
  )





condraw <- read.csv("C:/Users/mmahon/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Research/USGS Stream Macros/FISH/20201217.0745.FishSamp.csv",
                    colClasses = c("SiteNumber" = "character")) %>%
  mutate(SiteNumber = paste("USGS-",SiteNumber,sep = ""))

consub <- condraw %>% filter(!is.na(SpecificConductivity_µS)) %>%
  dplyr::select(SiteNumber, CollectionDate, SpecificConductivity_µS)

usgsfish <- usgsfish %>%
  left_join(consub %>% 
              mutate(SiteNumberCollectionDate = paste(SiteNumber, CollectionDate, sep = "_")) %>%
              dplyr::select(SiteNumberCollectionDate,SpecificConductivity_µS) %>%
              filter(SiteNumberCollectionDate %in% usgsfish$SiteNumberCollectionDate)
  )
usgsfish2 <- usgsfish2 %>%
  left_join(consub %>% 
              mutate(SiteNumberCollectionDate = paste(SiteNumber, CollectionDate, sep = "_")) %>%
              dplyr::select(SiteNumberCollectionDate,SpecificConductivity_µS) %>%
              filter(SiteNumberCollectionDate %in% usgsfish2$SiteNumberCollectionDate)
  )
usgsfish <- usgsfish %>%
  mutate(UseCond = ifelse(!is.na(SpecificConductivity_µS),
                          SpecificConductivity_µS,
                          Cond))
usgsfish2 <- usgsfish2 %>%
  mutate(UseCond = ifelse(!is.na(SpecificConductivity_µS),
                          SpecificConductivity_µS,
                          Cond))



##way to get up or downstream by buffer
# 
# final = data.frame()
# for(i in 180:nrow(missingsusgssites1)){
#   siteid <- missingsusgssites1$SiteNumber[i]
#   
#   nldiURLs <- list(UMwqp = paste("https://labs.waterdata.usgs.gov/api/nldi/linked-data/nwissite/",
#                                  siteid,
#                                  "/navigation/UM/wqp?distance=1",
#                                  sep = ""),
#                    DMwqp = paste("https://labs.waterdata.usgs.gov/api/nldi/linked-data/nwissite/",
#                                  siteid,
#                                  "/navigation/DM/wqp?distance=1",
#                                  sep = ""))
#   
#   nldi_data <- list()
#   for(n in names(nldiURLs)) {
#     nldi_data[n] <- list(sf::read_sf(nldiURLs[n][[1]]))
#     nldi_data$UMwqp$Location = "Upstream"
#     nldi_data$DMwqp$Location = "Downstream"
# 
#     hldr1 <- bind_rows(nldi_data) %>%
#       sf::st_drop_geometry() %>%
#       mutate(SITEID = siteid,
#              i = i)
#   }
#   final = bind_rows(final,hldr1)
# }

##skip errors: 12, 13, 27, 28, 86, 94, 132, 143, 144, 147, 148, 151, 153, 154,
##             161, 165, 168, 169, 178, 179, 


saveRDS(final,"nearbyUSGSsites.RDS")
final = readRDS("nearbyUSGSsites.RDS")

finalsub <- final %>%
  group_by(identifier, SITEID) %>%
  slice(1)

AllconductanceSecond <- TADA_DataRetrieval(characteristicName = c("Specific conductance",
                                                            "Specific conductivity",
                                                            "Conductivity",
                                                            "Conductance"),
                                     startDate = "1993-01-01",
                                     endDate = "2020-01-01",
                                     sampleMedia = c("Water", "water"),
                                     siteType = "Stream",
                                     siteid = finalsub$identifier,
                                     applyautoclean = F)


usgsfish.1 <- usgsfish %>%
  filter(is.na(UseCond)) %>%
  left_join(finalsub %>%
              dplyr::select(identifier, SITEID), by = join_by("SiteNumber" == "SITEID"))

usgsfish2.1 <- usgsfish2 %>%
  filter(is.na(UseCond)) %>%
  left_join(finalsub %>%
              dplyr::select(identifier, SITEID), by = join_by("SiteNumber" == "SITEID"))


allusgsfish.1 <- bind_rows(usgsfish.1 %>% dplyr::select(SiteNumber, identifier, CollectionDate, MonthEarlier, MonthLater),
                         usgsfish2.1 %>% dplyr::select(SiteNumber, identifier, CollectionDate, MonthEarlier, MonthLater)
) %>%
  group_by(SiteNumber,identifier, CollectionDate, MonthEarlier, MonthLater) %>%
  slice(1)



AllconductanceSecond2 = AllconductanceSecond %>%
  mutate(ActivityEndDate = as.Date(ActivityEndDate),
         ActivityStartDate = as.Date(ActivityStartDate)) %>%
  left_join(allusgsfish.1 %>% dplyr::select(identifier, CollectionDate, MonthEarlier, MonthLater),
            by = join_by("MonitoringLocationIdentifier" == "identifier"))

AllcondsubSecond <- AllconductanceSecond2 %>%
  filter(ActivityStartDate > MonthEarlier & ActivityStartDate <MonthLater |
           ActivityEndDate > MonthEarlier & ActivityEndDate <MonthLater) %>%
  dplyr::select(-MonthEarlier, -MonthLater)



AllcondsubSecond <- TADA_AutoClean(AllcondsubSecond)

AllcondsubSecondCleaned <- AllcondsubSecond %>%
  dplyr::select(ActivityStartDate,MonitoringLocationIdentifier,
                CollectionDate, TADA.ComparableDataIdentifier, 
                TADA.ResultMeasureValue, TADA.ResultMeasureValueDataTypes.Flag,
                ResultValueTypeName)

MeanConductivitySecond <- AllcondsubSecondCleaned %>%
  mutate(daydiff = as.numeric(ActivityStartDate - CollectionDate)) %>%
  group_by(MonitoringLocationIdentifier,CollectionDate) %>%
  mutate(minday = min(daydiff)) %>%
  ungroup() %>%
  filter(daydiff == minday) %>%
  group_by(MonitoringLocationIdentifier,CollectionDate) %>%
  summarize(Cond = mean(TADA.ResultMeasureValue, na.rm = T))

usgsfish.1 <- usgsfish.1 %>%
  mutate(identifierCollectionDate = paste(identifier, CollectionDate, sep = "_"))

usgsfish.1 <- usgsfish.1 %>%
  left_join(MeanConductivitySecond %>% 
              ungroup() %>%
              rename("SupCond" = "Cond") %>%
              mutate(SiteNumberCollectionDate = paste(MonitoringLocationIdentifier, CollectionDate, sep = "_")) %>%
              dplyr::select(-MonitoringLocationIdentifier, -CollectionDate) %>%
              filter(SiteNumberCollectionDate %in% usgsfish.1$identifierCollectionDate),
            by = join_by("identifierCollectionDate" == "SiteNumberCollectionDate")
  )

usgsfish2.1 <- usgsfish2.1 %>%
  mutate(identifierCollectionDate = paste(identifier, CollectionDate, sep = "_"))


usgsfish2.1 <- usgsfish2.1 %>%
  left_join(MeanConductivitySecond %>% 
              ungroup() %>%
              rename("SupCond" = "Cond") %>%
              mutate(SiteNumberCollectionDate = paste(MonitoringLocationIdentifier, CollectionDate, sep = "_")) %>%
              dplyr::select(-MonitoringLocationIdentifier, -CollectionDate) %>%
              filter(SiteNumberCollectionDate %in% usgsfish2.1$identifierCollectionDate),
            by = join_by("identifierCollectionDate" == "SiteNumberCollectionDate")
  )


usgsfish.1 <- usgsfish.1 %>%
  filter(!is.na(SupCond)) %>%
  group_by(SiteNumber, CollectionDate) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(SiteNumber,CollectionDate, SupCond)

usgsfish2.1 <- usgsfish2.1 %>%
  filter(!is.na(SupCond)) %>%
  group_by(SiteNumber, CollectionDate) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(SiteNumber,CollectionDate, SupCond)


usgsfish <- usgsfish %>%
  left_join(usgsfish.1 %>% 
              mutate(SiteNumberCollectionDate = paste(SiteNumber, CollectionDate, sep = "_")) %>%
              filter(SiteNumberCollectionDate %in% usgsfish$SiteNumberCollectionDate)
  )
usgsfish2 <- usgsfish2 %>%
  left_join(usgsfish2.1 %>% 
              mutate(SiteNumberCollectionDate = paste(SiteNumber, CollectionDate, sep = "_")) %>%
              filter(SiteNumberCollectionDate %in% usgsfish2$SiteNumberCollectionDate)
  )
usgsfish <- usgsfish %>%
  mutate(UseCond = ifelse(!is.na(UseCond),
                          UseCond,
                          SupCond))
usgsfish2 <- usgsfish2 %>%
  mutate(UseCond = ifelse(!is.na(UseCond),
                          UseCond,
                          SupCond))



sum(is.na(usgsfish2$UseCond)) / nrow(usgsfish2)
sum(is.na(usgsfish$UseCond)) / nrow(usgsfish)


saveRDS(usgsfish, "./Data/USGSOccur_Conductivity.rds")
saveRDS(usgsfish2, "./Data/USGSCPUE_Conductivity.rds")



##read in EPA conductivity dataset

fieldchem89 <- read.csv("https://www.epa.gov/sites/default/files/2015-11/fieldchemmeasure.csv")
labchem89 <- read.csv("https://www.epa.gov/sites/default/files/2015-09/chem.csv")

wholecond89 <- fieldchem89 %>%
  mutate(CONDUCTIVITY = as.numeric(gsub(" ", "", CONDUCTIVITY))) %>%
  dplyr::select(UID, SITE_ID, DATE_COL, CONDUCTIVITY) %>%
  left_join(labchem89 %>%
              dplyr::select(UID, SITE_ID, DATE_COL,COND,COND_ALERT)) %>%
  group_by(UID, SITE_ID, DATE_COL,COND_ALERT) %>%
  summarize(CONDUCTIVITY = mean(CONDUCTIVITY, na.rm = T),
            COND = mean(COND, na.rm = T))

wholecond89 <- wholecond89 %>%
  mutate(r1 = (COND / CONDUCTIVITY),
         r2 = ifelse((log10(r1) - floor(log10(r1))) >= 0.7,
                     floor(log10(r1)) + 1,
                     floor(log10(r1))),
         CONDUCTIVITY2 = CONDUCTIVITY * (10^r2),
         CONDUCTIVITY2 = ifelse(is.na(CONDUCTIVITY),
                                COND,
                                CONDUCTIVITY2))%>%
  mutate(DATE_COL = as.Date(DATE_COL,"%d-%b-%y"))

cor(wholecond89$COND,wholecond89$CONDUCTIVITY2,use = "complete.obs")


##1314

fieldchem1314 <- read.csv("https://www.epa.gov/sites/default/files/2019-04/nrsa1314_wide_field_meas_04232019.csv")
labchem1314 <- read.csv("https://www.epa.gov/sites/default/files/2019-04/nrsa1314_widechem_04232019.csv")

wholecond1314 <- fieldchem1314 %>%
  mutate(CONDUCTIVITY = as.numeric(gsub(" ", "", CONDUCTIVITY))) %>%
  dplyr::select(UID, SITE_ID, DATE_COL, CONDUCTIVITY,CONDUCTIVITY_FLAG) %>%
  left_join(labchem1314 %>%
              dplyr::select(UID, SITE_ID, DATE_COL,COND_RESULT,COND_NARS_FLAG,COND_QA_FLAG)) %>%
  group_by(UID, SITE_ID, DATE_COL,CONDUCTIVITY_FLAG,COND_NARS_FLAG,COND_QA_FLAG) %>%
  summarize(CONDUCTIVITY = mean(CONDUCTIVITY, na.rm = T),
            COND_RESULT = mean(COND_RESULT, na.rm = T))


wholecond1314 <- wholecond1314 %>%
  mutate(r1 = (COND_RESULT / CONDUCTIVITY),
         r2 = ifelse((log10(r1) - floor(log10(r1))) >= 0.7,
                     floor(log10(r1)) + 1,
                     floor(log10(r1))),
         CONDUCTIVITY2 = CONDUCTIVITY * (10^r2),
         CONDUCTIVITY2 = ifelse(is.na(CONDUCTIVITY),
                                COND_RESULT,
                                CONDUCTIVITY2)) %>%
  mutate(DATE_COL = as.Date(DATE_COL,"%m/%d/%Y"))
         
cor(wholecond1314$COND_RESULT,wholecond1314$CONDUCTIVITY2,use = "complete.obs")


  


##1819

fieldchem1819 <- read.csv("https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_field_wide.csv",fileEncoding="latin1")
colnames(fieldchem1819)[2] = "UID"
labchem1819 <- read.csv("https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_water_chemistry_chla_-_data.csv")

wholecond1819 <- fieldchem1819 %>%
  filter(SAMPLE_TYPE == "FIELDMEAS") %>%
  mutate(CONDUCTIVITY = as.numeric(gsub(" ", "", CONDUCTIVITY))) %>%
  dplyr::select(UID, SITE_ID, DATE_COL, CONDUCTIVITY) %>%
  left_join(labchem1819 %>%
              dplyr::select(UID, SITE_ID, DATE_COL,COND_RESULT,COND_NARS_FLAG,COND_QA_FLAG)) %>%
  group_by(UID, SITE_ID, DATE_COL,COND_NARS_FLAG,COND_QA_FLAG) %>%
  summarize(CONDUCTIVITY = mean(CONDUCTIVITY, na.rm = T),
            COND_RESULT = mean(COND_RESULT, na.rm = T))

wholecond1819 <- wholecond1819 %>%
  mutate(r1 = (COND_RESULT / CONDUCTIVITY),
         r2 = ifelse((log10(r1) - floor(log10(r1))) >= 0.7,
                     floor(log10(r1)) + 1,
                     floor(log10(r1))),
         CONDUCTIVITY2 = CONDUCTIVITY * (10^r2),
         CONDUCTIVITY2 = ifelse(is.na(CONDUCTIVITY),
                                COND_RESULT,
                                CONDUCTIVITY2)) %>%
  mutate(DATE_COL = as.Date(DATE_COL,"%m/%d/%Y"))

cor(wholecond1819$COND_RESULT,wholecond1819$CONDUCTIVITY2,use = "complete.obs")


EPA_wholecond <- bind_rows(wholecond89,wholecond1314,wholecond1819) %>%
  mutate(CONDUCTIVITY2 = ifelse(is.na(CONDUCTIVITY2),
                                CONDUCTIVITY,
                                CONDUCTIVITY2),
         Year = year(DATE_COL))

##save as RDS file
saveRDS(EPA_wholecond, "./Data/EPASites_Conductivity.rds")

fishdat1 <- fishdat1 %>%
  left_join(EPA_wholecond %>%
              ungroup() %>%
              group_by(UID) %>%
              slice(1) %>%
              ungroup() %>%
              mutate(UID = as.character(UID)) %>%
              dplyr::select(UID, CONDUCTIVITY2) %>%
              filter(UID %in% fishdat1$SampleID),
            by = join_by("SampleID" == "UID"))

fishdatCPUE <- fishdatCPUE %>%
  left_join(EPA_wholecond %>%
              ungroup() %>%
              group_by(UID) %>%
              slice(1) %>%
              ungroup() %>%
              mutate(UID = as.character(UID)) %>%
              dplyr::select(UID, CONDUCTIVITY2) %>%
              filter(UID %in% fishdat1$SampleID),
            by = join_by("SampleID" == "UID"))


fishdat1 = fishdat1 %>%
  mutate(CollectionDate = as.Date(CollectionDate)) %>%
  left_join(usgsfish %>%
              group_by(SiteNumber, CollectionDate) %>%
              slice(1) %>%
              ungroup() %>%
  dplyr::select(Agency:COMID, UseCond)) %>%
  mutate(WholeConductivity = ifelse(is.na(CONDUCTIVITY2),
                                    UseCond, 
                                    CONDUCTIVITY2))

fishdatCPUE = fishdatCPUE %>%
  mutate(CollectionDate = as.Date(CollectionDate)) %>%
  left_join(usgsfish2 %>%
              group_by(SiteNumber, CollectionDate) %>%
              slice(1) %>%
              ungroup() %>%
              dplyr::select(Agency:COMID, UseCond)) %>%
  mutate(WholeConductivity = ifelse(is.na(CONDUCTIVITY2),
                                    UseCond, 
                                    CONDUCTIVITY2))


write.csv(fishdatCPUE, "FishDivDatCPUEwCond.csv", row.names = F)
write.csv(fishdat1, "FishDivDat1.csv", row.names = F)

