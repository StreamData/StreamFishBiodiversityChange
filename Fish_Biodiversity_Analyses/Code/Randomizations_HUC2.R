##Randomizations

##within ecoregion
##need to calculate FDIS for each
library(vegan); library(fundiversity); library(svMisc)
##just need to break this into each individual ecoregion then; then rbind
##will need to break into separate dataframe first, then randomize, then 
##extract fdis(), then reconnect to ER dataframes, then rbind all ER dataframes
# nm1 <- nullmodel(fishdat[,27:738], "curveball")
# sm <- simulate(nm1, nsim=9999, seed = 32622)

#This then spits out 9999 iterations of the data; include observed as a fixed value
#then conduct fdis() on each one; go on from there
##Will likely need to remove those fish without traits before randomization!

##use deviation of fd_fdis function; accounts for lack of site names

fd_fdis <- function (traits, sp_com) 
{
  if (missing(traits) || is.null(traits)) {
    stop("Please provide a trait dataset", call. = FALSE)
  }
  if (is.data.frame(traits) || is.vector(traits)) {
    traits <- as.matrix(traits)
  }
  if (!is.numeric(traits)) {
    stop("Non-continuous trait data found in input traits. ", 
         "Please provide only continuous trait data", call. = FALSE)
  }
  traits <- fundiversity:::remove_species_without_trait(traits)
  if (!missing(sp_com)) {
    common_species <- fundiversity:::species_in_common(traits, sp_com)
    traits <- traits[common_species, , drop = FALSE]
    sp_com <- sp_com[, common_species, drop = FALSE]
  }
  else {
    sp_com <- matrix(1, ncol = nrow(traits), dimnames = list("s1", 
                                                             rownames(traits)))
  }
  site_abundances <- rowSums(sp_com, na.rm = TRUE)
  site_abundances[site_abundances == 0] <- 1
  sp_com <- sp_com/site_abundances
  centros <- sp_com %*% traits
  dists_centro <- future.apply::future_apply(centros, 1, function(centro) {
    sqrt(colSums((t(traits) - centro)^2))
  }, future.globals = FALSE)
  fdis_site <- diag(sp_com %*% dists_centro)
  
  data.frame(site = 1:nrow(sp_com), FDis = fdis_site, row.names = NULL)
}

##Whole community####

# fishdat1 <- fishdat1 %>%
#   left_join(flowdat %>% group_by(comid) %>%
#                                slice(1) %>% ungroup() %>%
#                                rename("COMID" = "comid"))
#split by Ecoregion
#need to generate 9 different datasets, remove all columns with colSums == 0
HUC01 = fishdat1 %>% filter(HUC2 == "01")
HUC01 = HUC01 %>%
  dplyr::select(-any_of(names(which(colSums(HUC01[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC01[,24:412]) != 0))))
HUC01 = HUC01[rowSums(HUC01[,46:ncol(HUC01)]) > 0,]

HUC02 = fishdat1 %>% filter(HUC2 == "02")
HUC02 = HUC02 %>%
  dplyr::select(-any_of(names(which(colSums(HUC02[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC02[,24:412]) != 0))))
HUC02 = HUC02[rowSums(HUC02[,46:ncol(HUC02)]) > 0,]

HUC03 = fishdat1 %>% filter(HUC2 == "03")
HUC03 = HUC03 %>%
  dplyr::select(-any_of(names(which(colSums(HUC03[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC03[,24:412]) != 0))))
HUC03 = HUC03[rowSums(HUC03[,46:ncol(HUC03)]) > 0,]

HUC04 = fishdat1 %>% filter(HUC2 == "04")
HUC04 = HUC04 %>%
  dplyr::select(-any_of(names(which(colSums(HUC04[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC04[,24:412]) != 0))))
HUC04 = HUC04[rowSums(HUC04[,46:ncol(HUC04)]) > 0,]

HUC05 = fishdat1 %>% filter(HUC2 == "05")
HUC05 = HUC05 %>%
  dplyr::select(-any_of(names(which(colSums(HUC05[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC05[,24:412]) != 0))))
HUC05 = HUC05[rowSums(HUC05[,46:ncol(HUC05)]) > 0,]

HUC06 = fishdat1 %>% filter(HUC2 == "06")
HUC06 = HUC06 %>%
  dplyr::select(-any_of(names(which(colSums(HUC06[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC06[,24:412]) != 0))))
HUC06 = HUC06[rowSums(HUC06[,46:ncol(HUC06)]) > 0,]

HUC07 = fishdat1 %>% filter(HUC2 == "07")
HUC07 = HUC07 %>%
  dplyr::select(-any_of(names(which(colSums(HUC07[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC07[,24:412]) != 0))))
HUC07 = HUC07[rowSums(HUC07[,46:ncol(HUC07)]) > 0,]

HUC08 = fishdat1 %>% filter(HUC2 == "08")
HUC08 = HUC08 %>%
  dplyr::select(-any_of(names(which(colSums(HUC08[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC08[,24:412]) != 0))))
HUC08 = HUC08[rowSums(HUC08[,46:ncol(HUC08)]) > 0,]

HUC09 = fishdat1 %>% filter(HUC2 == "09")
HUC09 = HUC09 %>%
  dplyr::select(-any_of(names(which(colSums(HUC09[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC09[,24:412]) != 0))))
HUC09 = HUC09[rowSums(HUC09[,46:ncol(HUC09)]) > 0,]

HUC10 = fishdat1 %>% filter(HUC2 == "10")
HUC10 = HUC10 %>%
  dplyr::select(-any_of(names(which(colSums(HUC10[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC10[,24:412]) != 0))))
HUC10 = HUC10[rowSums(HUC10[,46:ncol(HUC10)]) > 0,]

HUC11 = fishdat1 %>% filter(HUC2 == "11")
HUC11 = HUC11 %>%
  dplyr::select(-any_of(names(which(colSums(HUC11[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC11[,24:412]) != 0))))
HUC11 = HUC11[rowSums(HUC11[,46:ncol(HUC11)]) > 0,]

HUC12 = fishdat1 %>% filter(HUC2 == "12")
HUC12 = HUC12 %>%
  dplyr::select(-any_of(names(which(colSums(HUC12[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC12[,24:412]) != 0))))
HUC12 = HUC12[rowSums(HUC12[,46:ncol(HUC12)]) > 0,]

HUC13 = fishdat1 %>% filter(HUC2 == "13")
HUC13 = HUC13 %>%
  dplyr::select(-any_of(names(which(colSums(HUC13[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC13[,24:412]) != 0))))
HUC13 = HUC13[rowSums(HUC13[,46:ncol(HUC13)]) > 0,]

HUC14 = fishdat1 %>% filter(HUC2 == "14")
HUC14 = HUC14 %>%
  dplyr::select(-any_of(names(which(colSums(HUC14[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC14[,24:412]) != 0))))
HUC14 = HUC14[rowSums(HUC14[,46:ncol(HUC14)]) > 0,]

HUC15 = fishdat1 %>% filter(HUC2 == "15")
HUC15 = HUC15 %>%
  dplyr::select(-any_of(names(which(colSums(HUC15[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC15[,24:412]) != 0))))
HUC15 = HUC15[rowSums(HUC15[,46:ncol(HUC15)]) > 0,]

HUC16 = fishdat1 %>% filter(HUC2 == "16")
HUC16 = HUC16 %>%
  dplyr::select(-any_of(names(which(colSums(HUC16[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC16[,24:412]) != 0))))
HUC16 = HUC16[rowSums(HUC16[,46:ncol(HUC16)]) > 0,]

HUC17 = fishdat1 %>% filter(HUC2 == "17")
HUC17 = HUC17 %>%
  dplyr::select(-any_of(names(which(colSums(HUC17[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC17[,24:412]) != 0))))
HUC17 = HUC17[rowSums(HUC17[,46:ncol(HUC17)]) > 0,]

HUC18 = fishdat1 %>% filter(HUC2 == "18")
HUC18 = HUC18 %>%
  dplyr::select(-any_of(names(which(colSums(HUC18[,24:412]) == 0)))) %>%
  dplyr::select(Agency:MethodEffort, SppRichness:PredictedWettedWidth_km,
                any_of(names(which(colSums(HUC18[,24:412]) != 0))))
HUC18 = HUC18[rowSums(HUC18[,46:ncol(HUC18)]) > 0,]

##Before randonizations, need to filter these to only those spp present within
##the trait database, then filter out rows with no spp remaining




##What is happenning is that we are conducting separate burn-ins for 

##Gotelli and Ulrich recommendations:
##m (#row) x n (#col) initial burn-ins
##9999 simulations
##Indepenedent swaps OR 10 x n x m thins; thins are taking FOREVER,
##so go with 9999 independent swaps following n x m burn ins


nrow(HUC01)* ncol(HUC01[,-c(1:33)]) #
nrow(HUC02)* ncol(HUC02[,-c(1:33)]) #
nrow(HUC03)* ncol(HUC03[,-c(1:33)]) #
nrow(HUC04)* ncol(HUC04[,-c(1:33)]) #
nrow(HUC05)* ncol(HUC05[,-c(1:33)]) #
nrow(HUC06)* ncol(HUC06[,-c(1:33)]) #
nrow(HUC07)* ncol(HUC07[,-c(1:33)]) #
nrow(HUC08)* ncol(HUC08[,-c(1:33)]) #
nrow(HUC09)* ncol(HUC09[,-c(1:33)]) #
nrow(HUC10)* ncol(HUC10[,-c(1:33)]) #
nrow(HUC11)* ncol(HUC11[,-c(1:33)]) #
nrow(HUC12)* ncol(HUC12[,-c(1:33)]) #
nrow(HUC13)* ncol(HUC13[,-c(1:33)]) #
nrow(HUC14)* ncol(HUC14[,-c(1:33)]) #
nrow(HUC15)* ncol(HUC15[,-c(1:33)]) #
nrow(HUC16)* ncol(HUC16[,-c(1:33)]) #
nrow(HUC17)* ncol(HUC17[,-c(1:33)]) #
nrow(HUC18)* ncol(HUC18[,-c(1:33)]) #




##Generate 9999 different seeds, so that all iterations are independent;
##then run for loop to run simulate 9999 times

#01
set.seed(01)
HUC01seeds = sample(1:1000000,9999, replace = F)
HUC01hold = list()
for(i in 1:9999){
  HUC01hold[[i]] = simulate(nullmodel(HUC01[,-c(1:33)], "curveball"),
                          burnin = ifelse(nrow(HUC01)* ncol(HUC01[,-c(1:33)]) < 30000,
                                          30000,
                                          nrow(HUC01)* ncol(HUC01[,-c(1:33)]) ),
                                          
                          nsim = 1,
                          seed = HUC01seeds[i])
  progress(i,9999)
}
HUC01sm <- smbind(HUC01hold, MARGIN = 3)
saveRDS(HUC01sm, "HUC01sm.RDS")
HUC01sm = HUC01hold = NULL


#02
set.seed(02)
HUC02seeds = sample(1:1000000,9999, replace = F)
HUC02hold = list()
for(i in 1:9999){
  HUC02hold[[i]] = simulate(nullmodel(HUC02[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC02)* ncol(HUC02[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC02)* ncol(HUC02[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC02seeds[i])
  progress(i,9999)
}
HUC02sm <- smbind(HUC02hold, MARGIN = 3)
saveRDS(HUC02sm, "HUC02sm.RDS")
HUC02sm = HUC02hold = NULL

#03
set.seed(03)
HUC03seeds = sample(1:1000000,9999, replace = F)
HUC03hold = list()
for(i in 1:9999){
  HUC03hold[[i]] = simulate(nullmodel(HUC03[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC03)* ncol(HUC03[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC03)* ncol(HUC03[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC03seeds[i])
  progress(i,9999)
}
HUC03sm <- smbind(HUC03hold, MARGIN = 3)
saveRDS(HUC03sm, "HUC03sm.RDS")
HUC03sm = HUC03hold = NULL

#04
set.seed(04)
HUC04seeds = sample(1:1000000,9999, replace = F)
HUC04hold = list()
for(i in 1:9999){
  HUC04hold[[i]] = simulate(nullmodel(HUC04[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC04)* ncol(HUC04[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC04)* ncol(HUC04[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC04seeds[i])
  progress(i,9999)
}
HUC04sm <- smbind(HUC04hold, MARGIN = 3)
saveRDS(HUC04sm, "HUC04sm.RDS")
HUC04sm = HUC04hold = NULL

#05
set.seed(05)
HUC05seeds = sample(1:1000000,9999, replace = F)
HUC05hold = list()
for(i in 1:9999){
  HUC05hold[[i]] = simulate(nullmodel(HUC05[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC05)* ncol(HUC05[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC05)* ncol(HUC05[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC05seeds[i])
  progress(i,9999)
}
HUC05sm <- smbind(HUC05hold, MARGIN = 3)
saveRDS(HUC05sm, "HUC05sm.RDS")
HUC05sm = HUC05hold = NULL

#06
set.seed(06)
HUC06seeds = sample(1:1000000,9999, replace = F)
HUC06hold = list()
for(i in 1:9999){
  HUC06hold[[i]] = simulate(nullmodel(HUC06[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC06)* ncol(HUC06[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC06)* ncol(HUC06[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC06seeds[i])
  progress(i,9999)
}
HUC06sm <- smbind(HUC06hold, MARGIN = 3)
saveRDS(HUC06sm, "HUC06sm.RDS")
HUC06sm = HUC06hold = NULL

#07
set.seed(07)
HUC07seeds = sample(1:1000000,9999, replace = F)
HUC07hold = list()
for(i in 1:9999){
  HUC07hold[[i]] = simulate(nullmodel(HUC07[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC07)* ncol(HUC07[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC07)* ncol(HUC07[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC07seeds[i])
  progress(i,9999)
  
}
HUC07sm <- smbind(HUC07hold, MARGIN = 3)
saveRDS(HUC07sm, "HUC07sm.RDS")
HUC07sm = HUC07hold = NULL

#08
set.seed(08)
HUC08seeds = sample(1:1000000,9999, replace = F)
HUC08hold = list()
for(i in 1:9999){
  HUC08hold[[i]] = simulate(nullmodel(HUC08[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC08)* ncol(HUC08[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC08)* ncol(HUC08[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC08seeds[i])
  progress(i,9999)
}
HUC08sm <- smbind(HUC08hold, MARGIN = 3)
saveRDS(HUC08sm, "HUC08sm.RDS")
HUC08sm = HUC08hold = NULL

#09
set.seed(09)
HUC09seeds = sample(1:1000000,9999, replace = F)
HUC09hold = list()
for(i in 1:9999){
  HUC09hold[[i]] = simulate(nullmodel(HUC09[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC09)* ncol(HUC09[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC09)* ncol(HUC09[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC09seeds[i])
  progress(i,9999)
}
HUC09sm <- smbind(HUC09hold, MARGIN = 3)
saveRDS(HUC09sm, "HUC09sm.RDS")
HUC09sm = HUC09hold = NULL

#10
set.seed(10)
HUC10seeds = sample(1:1000000,9999, replace = F)
HUC10hold = list()
for(i in 1:9999){
  HUC10hold[[i]] = simulate(nullmodel(HUC10[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC10)* ncol(HUC10[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC10)* ncol(HUC10[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC10seeds[i])
  progress(i,9999)
}
HUC10sm <- smbind(HUC10hold, MARGIN = 3)
saveRDS(HUC10sm, "HUC10sm.RDS")
HUC10sm = HUC10hold = NULL

#11
set.seed(11)
HUC11seeds = sample(1:1000000,9999, replace = F)
HUC11hold = list()
for(i in 1:9999){
  HUC11hold[[i]] = simulate(nullmodel(HUC11[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC11)* ncol(HUC11[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC11)* ncol(HUC11[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC11seeds[i])
  progress(i,9999)
}
HUC11sm <- smbind(HUC11hold, MARGIN = 3)
saveRDS(HUC11sm, "HUC11sm.RDS")
HUC11sm = HUC11hold = NULL

#12
set.seed(12)
HUC12seeds = sample(1:1000000,9999, replace = F)
HUC12hold = list()
for(i in 1:9999){
  HUC12hold[[i]] = simulate(nullmodel(HUC12[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC12)* ncol(HUC12[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC12)* ncol(HUC12[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC12seeds[i])
  progress(i,9999)
}
HUC12sm <- smbind(HUC12hold, MARGIN = 3)
saveRDS(HUC12sm, "HUC12sm.RDS")
HUC12sm = HUC12hold = NULL

#13
set.seed(13)
HUC13seeds = sample(1:1000000,9999, replace = F)
HUC13hold = list()
for(i in 1:9999){
  HUC13hold[[i]] = simulate(nullmodel(HUC13[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC13)* ncol(HUC13[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC13)* ncol(HUC13[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC13seeds[i])
  progress(i,9999)
}
HUC13sm <- smbind(HUC13hold, MARGIN = 3)
saveRDS(HUC13sm, "HUC13sm.RDS")
HUC13sm = HUC13hold = NULL

#14
set.seed(14)
HUC14seeds = sample(1:1000000,9999, replace = F)
HUC14hold = list()
for(i in 1:9999){
  HUC14hold[[i]] = simulate(nullmodel(HUC14[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC14)* ncol(HUC14[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC14)* ncol(HUC14[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC14seeds[i])
  progress(i,9999)
}
HUC14sm <- smbind(HUC14hold, MARGIN = 3)
saveRDS(HUC14sm, "HUC14sm.RDS")
HUC14sm = HUC14hold = NULL

#15
set.seed(15)
HUC15seeds = sample(1:1000000,9999, replace = F)
HUC15hold = list()
for(i in 1:9999){
  HUC15hold[[i]] = simulate(nullmodel(HUC15[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC15)* ncol(HUC15[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC15)* ncol(HUC15[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC15seeds[i])
  progress(i,9999)
}
HUC15sm <- smbind(HUC15hold, MARGIN = 3)
saveRDS(HUC15sm, "HUC15sm.RDS")
HUC15sm = HUC15hold = NULL

#16
set.seed(16)
HUC16seeds = sample(1:1000000,9999, replace = F)
HUC16hold = list()
for(i in 1:9999){
  HUC16hold[[i]] = simulate(nullmodel(HUC16[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC16)* ncol(HUC16[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC16)* ncol(HUC16[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC16seeds[i])
  progress(i,9999)
}
HUC16sm <- smbind(HUC16hold, MARGIN = 3)
saveRDS(HUC16sm, "HUC16sm.RDS")
HUC16sm = HUC16hold = NULL

#17
set.seed(17)
HUC17seeds = sample(1:1000000,9999, replace = F)
HUC17hold = list()
for(i in 1:9999){
  HUC17hold[[i]] = simulate(nullmodel(HUC17[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC17)* ncol(HUC17[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC17)* ncol(HUC17[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC17seeds[i])
  progress(i,9999)
}
HUC17sm <- smbind(HUC17hold, MARGIN = 3)
saveRDS(HUC17sm, "HUC17sm.RDS")
HUC17sm = HUC17hold = NULL

#18
set.seed(18)
HUC18seeds = sample(1:1000000,9999, replace = F)
HUC18hold = list()
for(i in 1:9999){
  HUC18hold[[i]] = simulate(nullmodel(HUC18[,-c(1:33)], "curveball"),
                            burnin = ifelse(nrow(HUC18)* ncol(HUC18[,-c(1:33)]) < 30000,
                                            30000,
                                            nrow(HUC18)* ncol(HUC18[,-c(1:33)]) ),
                            
                            nsim = 1,
                            seed = HUC18seeds[i])
  progress(i,9999)
}
HUC18sm <- smbind(HUC18hold, MARGIN = 3)
saveRDS(HUC18sm, "HUC18sm.RDS")
HUC18sm = HUC18hold = NULL

#1,2,3 are done
#13,14,15,16,17,18 are done
#then run 6,8,9,12 DONE
#finally 4,5,7 DONE
# 11,10 NEED TO DO
##Next steps are to conduct fdis on each matrix; gather values, calculate means and SE

##Steps here:
#1) get trait dataset for ONLY those species in each region's species pool
#2) run gowdis on this; store it
#3) get species x site matrix only from each (should be easy to do, just get XXXsm[,,i],)
#3) A) but, will likely need to remove those taxa not in the trait database



##life histories first - start by redoing HUCs1-3; then change row names for the life history and habitat traits
# rownames(fulldatsub)[match(c("Lethenteron.appendix",
#                                     "Ericymba.amplamala",
#                                     "Etheostoma.chlorosomum"),
#                                   rownames(fulldatsub))] =
#   c("Lampetra.appendix",
#     "Notropis.amplamala",
#     "Etheostoma.chlorosoma")

##HUC01

HUC01sm <- readRDS("./Randomizations/HUC01sm.RDS")


HUC01traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC01[,-c(1:33)]),]

HUC01traitdis <- FD::gowdis(HUC01traits)
HUC01pcoa1 <- ape::pcoa(HUC01traitdis)
HUC01traitdispcoa <- as.matrix(HUC01pcoa1$vectors)


HUC01holdFD = matrix(nrow = nrow(HUC01), ncol = 9999)
for(i in 1:9999){
  HUC01holdFD[,i] = fd_fdis(traits = HUC01traitdispcoa, sp_com = as.matrix( HUC01sm[,colnames(HUC01sm[,,i]),i]))$FDis
  progress(i,9999)
}

HUC01fdis = fd_fdis(traits = HUC01traitdispcoa, sp_com = as.matrix(HUC01[,-c(1:33)]))$FDis
HUC01holdFD <- cbind(HUC01holdFD,c(HUC01fdis))
h1 <- (HUC01fdis - rowMeans(HUC01holdFD[,1:10000])) / apply(HUC01holdFD[,1:10000], 1, sd) 
HUC01$SES = h1
HUC01$fdis = HUC01fdis

saveRDS(HUC01, "./Randomizations/HUC01fdis.RDS")


##HUC02

HUC02sm <- readRDS("./Randomizations/HUC02sm.RDS")


HUC02traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC02[,-c(1:33)]),]

HUC02traitdis <- FD::gowdis(HUC02traits)
HUC02pcoa1 <- ape::pcoa(HUC02traitdis)
HUC02traitdispcoa <- as.matrix(HUC02pcoa1$vectors)


HUC02holdFD = matrix(nrow = nrow(HUC02), ncol = 9999)
for(i in 1:9999){
  HUC02holdFD[,i] = fd_fdis(traits = HUC02traitdispcoa, sp_com = as.matrix( HUC02sm[,colnames(HUC02sm[,,i]),i]))$FDis
  progress(i,9999)
}

HUC02fdis = fd_fdis(traits = HUC02traitdispcoa, sp_com = as.matrix(HUC02[,-c(1:33)]))$FDis
HUC02holdFD <- cbind(HUC02holdFD,c(HUC02fdis))
h1 <- (HUC02fdis - rowMeans(HUC02holdFD[,1:10000])) / apply(HUC02holdFD[,1:10000], 1, sd) 
HUC02$SES = h1
HUC02$fdis = HUC02fdis

saveRDS(HUC02, "./Randomizations/HUC02fdis.RDS")


##HUC03

HUC03sm <- readRDS("./Randomizations/HUC03sm.RDS")


HUC03traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC03[,-c(1:33)]),]



HUC03traitdis <- FD::gowdis(HUC03traits)
HUC03pcoa1 <- ape::pcoa(HUC03traitdis)
HUC03traitdispcoa <- as.matrix(HUC03pcoa1$vectors)


HUC03holdFD = matrix(nrow = nrow(HUC03), ncol = 9999)
for(i in 1:9999){
  HUC03holdFD[,i] = fd_fdis(traits = HUC03traitdispcoa, sp_com = as.matrix( HUC03sm[,colnames(HUC03sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC03fdis = fd_fdis(traits = HUC03traitdispcoa, sp_com = as.matrix(HUC03[,-c(1:33)]))$FDis
HUC03holdFD <- cbind(HUC03holdFD,c(HUC03fdis))
h1 <- (HUC03fdis - rowMeans(HUC03holdFD[,1:10000])) / apply(HUC03holdFD[,1:10000], 1, sd) 
HUC03$SES = h1
HUC03$fdis = HUC03fdis

saveRDS(HUC03, "./Randomizations/HUC03fdis.RDS")


##HUC04

HUC04sm <- readRDS("./Randomizations/HUC04sm.RDS")


HUC04traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC04[,-c(1:33)]),]

HUC04traitdis <- FD::gowdis(HUC04traits)
HUC04pcoa1 <- ape::pcoa(HUC04traitdis)
HUC04traitdispcoa <- as.matrix(HUC04pcoa1$vectors)


HUC04holdFD = matrix(nrow = nrow(HUC04), ncol = 9999)
for(i in 1:9999){
  HUC04holdFD[,i] = fd_fdis(traits = HUC04traitdispcoa, sp_com = as.matrix( HUC04sm[,colnames(HUC04sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC04fdis = fd_fdis(traits = HUC04traitdispcoa, sp_com = as.matrix(HUC04[,-c(1:33)]))$FDis
HUC04holdFD <- cbind(HUC04holdFD,c(HUC04fdis))
h1 <- (HUC04fdis - rowMeans(HUC04holdFD[,1:10000])) / apply(HUC04holdFD[,1:10000], 1, sd) 
HUC04$SES = h1
HUC04$fdis = HUC04fdis

saveRDS(HUC04, "./Randomizations/HUC04fdis.RDS")


##HUC05

HUC05sm <- readRDS("./Randomizations/HUC05sm.RDS")


HUC05traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC05[,-c(1:33)]),]

HUC05traitdis <- FD::gowdis(HUC05traits)
HUC05pcoa1 <- ape::pcoa(HUC05traitdis)
HUC05traitdispcoa <- as.matrix(HUC05pcoa1$vectors)


HUC05holdFD = matrix(nrow = nrow(HUC05), ncol = 9999)
for(i in 1:9999){
  HUC05holdFD[,i] = fd_fdis(traits = HUC05traitdispcoa, sp_com = as.matrix( HUC05sm[,colnames(HUC05sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC05fdis = fd_fdis(traits = HUC05traitdispcoa, sp_com = as.matrix(HUC05[,-c(1:33)]))$FDis
HUC05holdFD <- cbind(HUC05holdFD,c(HUC05fdis))
h1 <- (HUC05fdis - rowMeans(HUC05holdFD[,1:10000])) / apply(HUC05holdFD[,1:10000], 1, sd) 
HUC05$SES = h1
HUC05$fdis = HUC05fdis

saveRDS(HUC05, "./Randomizations/HUC05fdis.RDS")

##HUC06

HUC06sm <- readRDS("./Randomizations/HUC06sm.RDS")


HUC06traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC06[,-c(1:33)]),]

HUC06traitdis <- FD::gowdis(HUC06traits)
HUC06pcoa1 <- ape::pcoa(HUC06traitdis)
HUC06traitdispcoa <- as.matrix(HUC06pcoa1$vectors)


HUC06holdFD = matrix(nrow = nrow(HUC06), ncol = 9999)
for(i in 1:9999){
  HUC06holdFD[,i] = fd_fdis(traits = HUC06traitdispcoa, sp_com = as.matrix( HUC06sm[,colnames(HUC06sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC06fdis = fd_fdis(traits = HUC06traitdispcoa, sp_com = as.matrix(HUC06[,-c(1:33)]))$FDis
HUC06holdFD <- cbind(HUC06holdFD,c(HUC06fdis))
h1 <- (HUC06fdis - rowMeans(HUC06holdFD[,1:10000])) / apply(HUC06holdFD[,1:10000], 1, sd) 
HUC06$SES = h1
HUC06$fdis = HUC06fdis

saveRDS(HUC06, "./Randomizations/HUC06fdis.RDS")


##HUC07

HUC07sm <- readRDS("./Randomizations/HUC07sm.RDS")


HUC07traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC07[,-c(1:33)]),]

HUC07traitdis <- FD::gowdis(HUC07traits)
HUC07pcoa1 <- ape::pcoa(HUC07traitdis)
HUC07traitdispcoa <- as.matrix(HUC07pcoa1$vectors)


HUC07holdFD = matrix(nrow = nrow(HUC07), ncol = 9999)
for(i in 1:9999){
  HUC07holdFD[,i] = fd_fdis(traits = HUC07traitdispcoa, sp_com = as.matrix( HUC07sm[,colnames(HUC07sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC07fdis = fd_fdis(traits = HUC07traitdispcoa, sp_com = as.matrix(HUC07[,-c(1:33)]))$FDis
HUC07holdFD <- cbind(HUC07holdFD,c(HUC07fdis))
h1 <- (HUC07fdis - rowMeans(HUC07holdFD[,1:10000])) / apply(HUC07holdFD[,1:10000], 1, sd) 
HUC07$SES = h1
HUC07$fdis = HUC07fdis

saveRDS(HUC07, "./Randomizations/HUC07fdis.RDS")


##HUC08

HUC08sm <- readRDS("./Randomizations/HUC08sm.RDS")


HUC08traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC08[,-c(1:33)]),]

HUC08traitdis <- FD::gowdis(HUC08traits)
HUC08pcoa1 <- ape::pcoa(HUC08traitdis)
HUC08traitdispcoa <- as.matrix(HUC08pcoa1$vectors)


HUC08holdFD = matrix(nrow = nrow(HUC08), ncol = 9999)
for(i in 1:9999){
  HUC08holdFD[,i] = fd_fdis(traits = HUC08traitdispcoa, sp_com = as.matrix( HUC08sm[,colnames(HUC08sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC08fdis = fd_fdis(traits = HUC08traitdispcoa, sp_com = as.matrix(HUC08[,-c(1:33)]))$FDis
HUC08holdFD <- cbind(HUC08holdFD,c(HUC08fdis))
h1 <- (HUC08fdis - rowMeans(HUC08holdFD[,1:10000])) / apply(HUC08holdFD[,1:10000], 1, sd) 
HUC08$SES = h1
HUC08$fdis = HUC08fdis

saveRDS(HUC08, "./Randomizations/HUC08fdis.RDS")


##HUC09

HUC09sm <- readRDS("./Randomizations/HUC09sm.RDS")


HUC09traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC09[,-c(1:33)]),]

HUC09traitdis <- FD::gowdis(HUC09traits)
HUC09pcoa1 <- ape::pcoa(HUC09traitdis)
HUC09traitdispcoa <- as.matrix(HUC09pcoa1$vectors)


HUC09holdFD = matrix(nrow = nrow(HUC09), ncol = 9999)
for(i in 1:9999){
  HUC09holdFD[,i] = fd_fdis(traits = HUC09traitdispcoa, sp_com = as.matrix( HUC09sm[,colnames(HUC09sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC09fdis = fd_fdis(traits = HUC09traitdispcoa, sp_com = as.matrix(HUC09[,-c(1:33)]))$FDis
HUC09holdFD <- cbind(HUC09holdFD,c(HUC09fdis))
h1 <- (HUC09fdis - rowMeans(HUC09holdFD[,1:10000])) / apply(HUC09holdFD[,1:10000], 1, sd) 
HUC09$SES = h1
HUC09$fdis = HUC09fdis

saveRDS(HUC09, "./Randomizations/HUC09fdis.RDS")


##HUC10

HUC10sm <- readRDS("./Randomizations/HUC10sm.RDS")


HUC10traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC10[,-c(1:33)]),]

HUC10traitdis <- FD::gowdis(HUC10traits)
HUC10pcoa1 <- ape::pcoa(HUC10traitdis)
HUC10traitdispcoa <- as.matrix(HUC10pcoa1$vectors)


HUC10holdFD = matrix(nrow = nrow(HUC10), ncol = 9999)
for(i in 1:9999){
  HUC10holdFD[,i] = fd_fdis(traits = HUC10traitdispcoa, sp_com = as.matrix( HUC10sm[,colnames(HUC10sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC10fdis = fd_fdis(traits = HUC10traitdispcoa, sp_com = as.matrix(HUC10[,-c(1:33)]))$FDis
HUC10holdFD <- cbind(HUC10holdFD,c(HUC10fdis))
h1 <- (HUC10fdis - rowMeans(HUC10holdFD[,1:10000])) / apply(HUC10holdFD[,1:10000], 1, sd) 
HUC10$SES = h1
HUC10$fdis = HUC10fdis


saveRDS(HUC10, "./Randomizations/HUC10fdis.RDS")

##HUC11

HUC11sm <- readRDS("./Randomizations/HUC11sm.RDS")


HUC11traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC11[,-c(1:33)]),]

HUC11traitdis <- FD::gowdis(HUC11traits)
HUC11pcoa1 <- ape::pcoa(HUC11traitdis)
HUC11traitdispcoa <- as.matrix(HUC11pcoa1$vectors)


HUC11holdFD = matrix(nrow = nrow(HUC11), ncol = 9999)
for(i in 1:9999){
  HUC11holdFD[,i] = fd_fdis(traits = HUC11traitdispcoa, sp_com = as.matrix( HUC11sm[,colnames(HUC11sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC11fdis = fd_fdis(traits = HUC11traitdispcoa, sp_com = as.matrix(HUC11[,-c(1:33)]))$FDis
HUC11holdFD <- cbind(HUC11holdFD,c(HUC11fdis))
h1 <- (HUC11fdis - rowMeans(HUC11holdFD[,1:10000])) / apply(HUC11holdFD[,1:10000], 1, sd) 
HUC11$SES = h1
HUC11$fdis = HUC11fdis

saveRDS(HUC11, "./Randomizations/HUC11fdis.RDS")


##HUC12

HUC12sm <- readRDS("./Randomizations/HUC12sm.RDS")


HUC12traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC12[,-c(1:33)]),]

HUC12traitdis <- FD::gowdis(HUC12traits)
HUC12pcoa1 <- ape::pcoa(HUC12traitdis)
HUC12traitdispcoa <- as.matrix(HUC12pcoa1$vectors)


HUC12holdFD = matrix(nrow = nrow(HUC12), ncol = 9999)
for(i in 1:9999){
  HUC12holdFD[,i] = fd_fdis(traits = HUC12traitdispcoa, sp_com = as.matrix( HUC12sm[,colnames(HUC12sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC12fdis = fd_fdis(traits = HUC12traitdispcoa, sp_com = as.matrix(HUC12[,-c(1:33)]))$FDis
HUC12holdFD <- cbind(HUC12holdFD,c(HUC12fdis))
h1 <- (HUC12fdis - rowMeans(HUC12holdFD[,1:10000])) / apply(HUC12holdFD[,1:10000], 1, sd) 
HUC12$SES = h1
HUC12$fdis = HUC12fdis

saveRDS(HUC12, "./Randomizations/HUC12fdis.RDS")


##HUC13

HUC13sm <- readRDS("./Randomizations/HUC13sm.RDS")


HUC13traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC13[,-c(1:33)]),]

HUC13traitdis <- FD::gowdis(HUC13traits)
HUC13pcoa1 <- ape::pcoa(HUC13traitdis)
HUC13traitdispcoa <- as.matrix(HUC13pcoa1$vectors)


HUC13holdFD = matrix(nrow = nrow(HUC13), ncol = 9999)
for(i in 1:9999){
  HUC13holdFD[,i] = fd_fdis(traits = HUC13traitdispcoa, sp_com = as.matrix( HUC13sm[,colnames(HUC13sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC13fdis = fd_fdis(traits = HUC13traitdispcoa, sp_com = as.matrix(HUC13[,-c(1:33)]))$FDis
HUC13holdFD <- cbind(HUC13holdFD,c(HUC13fdis))
h1 <- (HUC13fdis - rowMeans(HUC13holdFD[,1:10000])) / apply(HUC13holdFD[,1:10000], 1, sd) 
HUC13$SES = h1
HUC13$fdis = HUC13fdis

saveRDS(HUC13, "./Randomizations/HUC13fdis.RDS")


##HUC14

HUC14sm <- readRDS("./Randomizations/HUC14sm.RDS")


HUC14traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC14[,-c(1:33)]),]

HUC14traitdis <- FD::gowdis(HUC14traits)
HUC14pcoa1 <- ape::pcoa(HUC14traitdis)
HUC14traitdispcoa <- as.matrix(HUC14pcoa1$vectors)


HUC14holdFD = matrix(nrow = nrow(HUC14), ncol = 9999)
for(i in 1:9999){
  HUC14holdFD[,i] = fd_fdis(traits = HUC14traitdispcoa, sp_com = as.matrix( HUC14sm[,colnames(HUC14sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC14fdis = fd_fdis(traits = HUC14traitdispcoa, sp_com = as.matrix(HUC14[,-c(1:33)]))$FDis
HUC14holdFD <- cbind(HUC14holdFD,c(HUC14fdis))
h1 <- (HUC14fdis - rowMeans(HUC14holdFD[,1:10000])) / apply(HUC14holdFD[,1:10000], 1, sd) 
HUC14$SES = h1
HUC14$fdis = HUC14fdis

saveRDS(HUC14, "./Randomizations/HUC14fdis.RDS")


##HUC15

HUC15sm <- readRDS("./Randomizations/HUC15sm.RDS")


HUC15traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC15[,-c(1:33)]),]

HUC15traitdis <- FD::gowdis(HUC15traits)
HUC15pcoa1 <- ape::pcoa(HUC15traitdis)
HUC15traitdispcoa <- as.matrix(HUC15pcoa1$vectors)


HUC15holdFD = matrix(nrow = nrow(HUC15), ncol = 9999)
for(i in 1:9999){
  HUC15holdFD[,i] = fd_fdis(traits = HUC15traitdispcoa, sp_com = as.matrix( HUC15sm[,colnames(HUC15sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC15fdis = fd_fdis(traits = HUC15traitdispcoa, sp_com = as.matrix(HUC15[,-c(1:33)]))$FDis
HUC15holdFD <- cbind(HUC15holdFD,c(HUC15fdis))
h1 <- (HUC15fdis - rowMeans(HUC15holdFD[,1:10000])) / apply(HUC15holdFD[,1:10000], 1, sd) 
HUC15$SES = h1
HUC15$fdis = HUC15fdis

saveRDS(HUC15, "./Randomizations/HUC15fdis.RDS")


##HUC16

HUC16sm <- readRDS("./Randomizations/HUC16sm.RDS")


HUC16traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC16[,-c(1:33)]),]

HUC16traitdis <- FD::gowdis(HUC16traits)
HUC16pcoa1 <- ape::pcoa(HUC16traitdis)
HUC16traitdispcoa <- as.matrix(HUC16pcoa1$vectors)


HUC16holdFD = matrix(nrow = nrow(HUC16), ncol = 9999)
for(i in 1:9999){
  HUC16holdFD[,i] = fd_fdis(traits = HUC16traitdispcoa, sp_com = as.matrix( HUC16sm[,colnames(HUC16sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC16fdis = fd_fdis(traits = HUC16traitdispcoa, sp_com = as.matrix(HUC16[,-c(1:33)]))$FDis
HUC16holdFD <- cbind(HUC16holdFD,c(HUC16fdis))
h1 <- (HUC16fdis - rowMeans(HUC16holdFD[,1:10000])) / apply(HUC16holdFD[,1:10000], 1, sd) 
HUC16$SES = h1
HUC16$fdis = HUC16fdis

saveRDS(HUC16, "./Randomizations/HUC16fdis.RDS")


##HUC17

HUC17sm <- readRDS("./Randomizations/HUC17sm.RDS")


HUC17traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC17[,-c(1:33)]),]

HUC17traitdis <- FD::gowdis(HUC17traits)
HUC17pcoa1 <- ape::pcoa(HUC17traitdis)
HUC17traitdispcoa <- as.matrix(HUC17pcoa1$vectors)


HUC17holdFD = matrix(nrow = nrow(HUC17), ncol = 9999)
for(i in 1:9999){
  HUC17holdFD[,i] = fd_fdis(traits = HUC17traitdispcoa, sp_com = as.matrix( HUC17sm[,colnames(HUC17sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC17fdis = fd_fdis(traits = HUC17traitdispcoa, sp_com = as.matrix(HUC17[,-c(1:33)]))$FDis
HUC17holdFD <- cbind(HUC17holdFD,c(HUC17fdis))
h1 <- (HUC17fdis - rowMeans(HUC17holdFD[,1:10000])) / apply(HUC17holdFD[,1:10000], 1, sd) 
HUC17$SES = h1
HUC17$fdis = HUC17fdis

saveRDS(HUC17, "./Randomizations/HUC17fdis.RDS")


##HUC18

HUC18sm <- readRDS("./Randomizations/HUC18sm.RDS")


HUC18traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC18[,-c(1:33)]),]

HUC18traitdis <- FD::gowdis(HUC18traits)
HUC18pcoa1 <- ape::pcoa(HUC18traitdis)
HUC18traitdispcoa <- as.matrix(HUC18pcoa1$vectors)


HUC18holdFD = matrix(nrow = nrow(HUC18), ncol = 9999)
for(i in 1:9999){
  HUC18holdFD[,i] = fd_fdis(traits = HUC18traitdispcoa, sp_com = as.matrix( HUC18sm[,colnames(HUC18sm[,,i]),i]))$FDis
  progress(i, 9999)
}

HUC18fdis = fd_fdis(traits = HUC18traitdispcoa, sp_com = as.matrix(HUC18[,-c(1:33)]))$FDis
HUC18holdFD <- cbind(HUC18holdFD,c(HUC18fdis))
h1 <- (HUC18fdis - rowMeans(HUC18holdFD[,1:10000])) / apply(HUC18holdFD[,1:10000], 1, sd) 
HUC18$SES = h1
HUC18$fdis = HUC18fdis

saveRDS(HUC18, "./Randomizations/HUC18fdis.RDS")


##the following need to be done:
# 7, 10


HUC01 <- readRDS("./Randomizations/HUC01fdis.RDS");
HUC02 <- readRDS("./Randomizations/HUC02fdis.RDS")
HUC03 <- readRDS("./Randomizations/HUC03fdis.RDS");HUC04 <- readRDS("./Randomizations/HUC04fdis.RDS")
HUC05 <- readRDS("./Randomizations/HUC05fdis.RDS");HUC06 <- readRDS("./Randomizations/HUC06fdis.RDS")
HUC07 <- readRDS("./Randomizations/HUC07fdis.RDS");HUC08 <- readRDS("./Randomizations/HUC08fdis.RDS")
HUC09 <- readRDS("./Randomizations/HUC09fdis.RDS");HUC10 <- readRDS("./Randomizations/HUC10fdis.RDS")
HUC11 <- readRDS("./Randomizations/HUC11fdis.RDS");HUC12 <- readRDS("./Randomizations/HUC12fdis.RDS")
HUC13 <- readRDS("./Randomizations/HUC13fdis.RDS");HUC14 <- readRDS("./Randomizations/HUC14fdis.RDS")
HUC15 <- readRDS("./Randomizations/HUC15fdis.RDS");HUC16 <- readRDS("./Randomizations/HUC16fdis.RDS")
HUC17 <- readRDS("./Randomizations/HUC17fdis.RDS");HUC18 <- readRDS("./Randomizations/HUC18fdis.RDS")

WHOLEHUC = bind_rows(HUC01[,c(1:33,ncol(HUC01)-1,ncol(HUC01))] %>% mutate(HUC2 = "01"),
                     HUC02[,c(1:33,ncol(HUC02)-1,ncol(HUC02))]%>% mutate(HUC2 = "02"),
                   HUC03[,c(1:33,ncol(HUC03)-1,ncol(HUC03))]%>% mutate(HUC2 = "03"),
                   HUC04[,c(1:33,ncol(HUC04)-1,ncol(HUC04))]%>% mutate(HUC2 = "04"),
                   HUC05[,c(1:33,ncol(HUC05)-1,ncol(HUC05))]%>% mutate(HUC2 = "05"),
                   HUC06[,c(1:33,ncol(HUC06)-1,ncol(HUC06))]%>% mutate(HUC2 = "06"),
                   HUC07[,c(1:33,ncol(HUC07)-1,ncol(HUC07))]%>% mutate(HUC2 = "07"),
                   HUC08[,c(1:33,ncol(HUC08)-1,ncol(HUC08))]%>% mutate(HUC2 = "08"),
                   HUC09[,c(1:33,ncol(HUC09)-1,ncol(HUC09))]%>% mutate(HUC2 = "09"), 
                   HUC10[,c(1:33,ncol(HUC10)-1,ncol(HUC10))]%>% mutate(HUC2 = "10"),
                   HUC11[,c(1:33,ncol(HUC11)-1,ncol(HUC11))]%>% mutate(HUC2 = "11"),
                   HUC12[,c(1:33,ncol(HUC12)-1,ncol(HUC12))]%>% mutate(HUC2 = "12"),
                   HUC13[,c(1:33,ncol(HUC13)-1,ncol(HUC13))]%>% mutate(HUC2 = "13"),
                   HUC14[,c(1:33,ncol(HUC14)-1,ncol(HUC14))]%>% mutate(HUC2 = "14"),
                   HUC15[,c(1:33,ncol(HUC15)-1,ncol(HUC15))]%>% mutate(HUC2 = "15"), 
                   HUC16[,c(1:33,ncol(HUC16)-1,ncol(HUC16))]%>% mutate(HUC2 = "16"),
                   HUC17[,c(1:33,ncol(HUC17)-1,ncol(HUC17))]%>% mutate(HUC2 = "17"), 
                   HUC18[,c(1:33,ncol(HUC18)-1,ncol(HUC18))]%>% mutate(HUC2 = "18")
                   )

saveRDS(WHOLEHUC,"SESWholeCommFinal.RDS")

##native####

##will need to filter species to only non-managed ones - need to do this in
##multiple spots for each HUC
gamey = read.csv("./Data/GameFishDesignation.csv")
nonnative1 = read.csv("./Data/USGSinput2_update_screened_nonnative.csv")
nonnative2 = read.csv("./Data/USGSinput1_update_screened_nonnative.csv")
nonnative = rbind(nonnative1, nonnative2)

nnh2 = nonnative %>%
  mutate(HUC2 = substr(str_pad(as.character(HUC8), 8, pad = "0"),1,2)) %>%
  dplyr::select(Scientific.Name,HUC2, Native) %>%
  filter(Native == "No") %>%
  mutate(Scientific.Name = gsub(" ","\\.",Scientific.Name)) %>%
  group_by(Scientific.Name, HUC2) %>%
  slice(1) %>%
  ungroup()

nnh2 = split(nnh2, ~HUC2)

##1,2,4,6,8,9, 13-18, 5,7,12 11,10 are done
## 3 running

##HUC01

HUC01sm <- readRDS("./Randomizations/HUC01sm.RDS")

HUC01traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC01[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["01"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                      "SES", "fdis")))),]

HUC01traitdis <- FD::gowdis(HUC01traits)
HUC01pcoa1 <- ape::pcoa(HUC01traitdis)
HUC01traitdispcoa <- as.matrix(HUC01pcoa1$vectors)

HUC01holdFD = matrix(nrow = nrow(HUC01), ncol = 9999)
for(i in 1:9999){
  
  HUC01holdFD[,i] = fd_fdis(traits = HUC01traitdispcoa,
                            sp_com = as.matrix(HUC01sm[,!(colnames(HUC01sm[,,i]) %in% unique(c(nnh2[["01"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC01fdis = fd_fdis(traits = HUC01traitdispcoa, sp_com = as.matrix(HUC01[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["01"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC01holdFD <- cbind(HUC01holdFD,c(HUC01fdis))
h1 <- (HUC01fdis - rowMeans(HUC01holdFD[,1:10000])) / apply(HUC01holdFD[,1:10000], 1, sd) 
HUC01$SES = h1
HUC01$fdis = HUC01fdis

saveRDS(HUC01, "./Randomizations/HUC01fdis_nng.RDS")
HUC01sm = NULL

##HUC02


HUC02sm <- readRDS("./Randomizations/HUC02sm.RDS")

HUC02traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC02[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["02"]]$Scientific.Name, gsub(" ","\\.",gamey$Species, "Cyprinus.carpio"))),
                                                                                        "SES", "fdis")))),]

HUC02traitdis <- FD::gowdis(HUC02traits)
HUC02pcoa1 <- ape::pcoa(HUC02traitdis)
HUC02traitdispcoa <- as.matrix(HUC02pcoa1$vectors)

HUC02holdFD = matrix(nrow = nrow(HUC02), ncol = 9999)
for(i in 1:9999){
  
  HUC02holdFD[,i] = fd_fdis(traits = HUC02traitdispcoa,
                            sp_com = as.matrix(HUC02sm[,!(colnames(HUC02sm[,,i]) %in% unique(c(nnh2[["02"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC02fdis = fd_fdis(traits = HUC02traitdispcoa, sp_com = as.matrix(HUC02[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["02"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC02holdFD <- cbind(HUC02holdFD,c(HUC02fdis))
h1 <- (HUC02fdis - rowMeans(HUC02holdFD[,1:10000])) / apply(HUC02holdFD[,1:10000], 1, sd) 
HUC02$SES = h1
HUC02$fdis = HUC02fdis

saveRDS(HUC02, "./Randomizations/HUC02fdis_nng.RDS")
HUC02sm = NULL


##HUC03
HUC03sm <- readRDS("./Randomizations/HUC03sm.RDS")

HUC03traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC03[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["03"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC03traitdis <- FD::gowdis(HUC03traits)
HUC03pcoa1 <- ape::pcoa(HUC03traitdis)
HUC03traitdispcoa <- as.matrix(HUC03pcoa1$vectors)

HUC03holdFD = matrix(nrow = nrow(HUC03), ncol = 9999)
for(i in 1:9999){
  
  HUC03holdFD[,i] = fd_fdis(traits = HUC03traitdispcoa,
                            sp_com = as.matrix(HUC03sm[,!(colnames(HUC03sm[,,i]) %in% unique(c(nnh2[["03"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC03fdis = fd_fdis(traits = HUC03traitdispcoa, sp_com = as.matrix(HUC03[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["03"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC03holdFD <- cbind(HUC03holdFD,c(HUC03fdis))
h1 <- (HUC03fdis - rowMeans(HUC03holdFD[,1:10000])) / apply(HUC03holdFD[,1:10000], 1, sd) 
HUC03$SES = h1
HUC03$fdis = HUC03fdis

saveRDS(HUC03, "./Randomizations/HUC03fdis_nng.RDS")
HUC03sm = NULL


##HUC04
HUC04sm <- readRDS("./Randomizations/HUC04sm.RDS")

HUC04traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC04[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["04"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC04traitdis <- FD::gowdis(HUC04traits)
HUC04pcoa1 <- ape::pcoa(HUC04traitdis)
HUC04traitdispcoa <- as.matrix(HUC04pcoa1$vectors)

HUC04holdFD = matrix(nrow = nrow(HUC04), ncol = 9999)
for(i in 1:9999){
  
  HUC04holdFD[,i] = fd_fdis(traits = HUC04traitdispcoa,
                            sp_com = as.matrix(HUC04sm[,!(colnames(HUC04sm[,,i]) %in% unique(c(nnh2[["04"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC04fdis = fd_fdis(traits = HUC04traitdispcoa, sp_com = as.matrix(HUC04[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["04"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC04holdFD <- cbind(HUC04holdFD,c(HUC04fdis))
h1 <- (HUC04fdis - rowMeans(HUC04holdFD[,1:10000])) / apply(HUC04holdFD[,1:10000], 1, sd) 
HUC04$SES = h1
HUC04$fdis = HUC04fdis

saveRDS(HUC04, "./Randomizations/HUC04fdis_nng.RDS")
HUC04sm = NULL


##HUC05
HUC05sm <- readRDS("./Randomizations/HUC05sm.RDS")

HUC05traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC05[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["05"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC05traitdis <- FD::gowdis(HUC05traits)
HUC05pcoa1 <- ape::pcoa(HUC05traitdis)
HUC05traitdispcoa <- as.matrix(HUC05pcoa1$vectors)

HUC05holdFD = matrix(nrow = nrow(HUC05), ncol = 9999)
for(i in 1:9999){
  
  HUC05holdFD[,i] = fd_fdis(traits = HUC05traitdispcoa,
                            sp_com = as.matrix(HUC05sm[,!(colnames(HUC05sm[,,i]) %in% unique(c(nnh2[["05"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC05fdis = fd_fdis(traits = HUC05traitdispcoa, sp_com = as.matrix(HUC05[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["05"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC05holdFD <- cbind(HUC05holdFD,c(HUC05fdis))
h1 <- (HUC05fdis - rowMeans(HUC05holdFD[,1:10000])) / apply(HUC05holdFD[,1:10000], 1, sd) 
HUC05$SES = h1
HUC05$fdis = HUC05fdis

saveRDS(HUC05, "./Randomizations/HUC05fdis_nng.RDS")
HUC05sm = NULL

##HUC06
HUC06sm <- readRDS("./Randomizations/HUC06sm.RDS")

HUC06traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC06[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["06"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC06traitdis <- FD::gowdis(HUC06traits)
HUC06pcoa1 <- ape::pcoa(HUC06traitdis)
HUC06traitdispcoa <- as.matrix(HUC06pcoa1$vectors)

HUC06holdFD = matrix(nrow = nrow(HUC06), ncol = 9999)
for(i in 1:9999){
  
  HUC06holdFD[,i] = fd_fdis(traits = HUC06traitdispcoa,
                            sp_com = as.matrix(HUC06sm[,!(colnames(HUC06sm[,,i]) %in% unique(c(nnh2[["06"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC06fdis = fd_fdis(traits = HUC06traitdispcoa, sp_com = as.matrix(HUC06[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["06"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC06holdFD <- cbind(HUC06holdFD,c(HUC06fdis))
h1 <- (HUC06fdis - rowMeans(HUC06holdFD[,1:10000])) / apply(HUC06holdFD[,1:10000], 1, sd) 
HUC06$SES = h1
HUC06$fdis = HUC06fdis

saveRDS(HUC06, "./Randomizations/HUC06fdis_nng.RDS")
HUC06sm = NULL


##HUC07
HUC07sm <- readRDS("./Randomizations/HUC07sm.RDS")

HUC07traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC07[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["07"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC07traitdis <- FD::gowdis(HUC07traits)
HUC07pcoa1 <- ape::pcoa(HUC07traitdis)
HUC07traitdispcoa <- as.matrix(HUC07pcoa1$vectors)

HUC07holdFD = matrix(nrow = nrow(HUC07), ncol = 9999)
for(i in 1:9999){
  
  HUC07holdFD[,i] = fd_fdis(traits = HUC07traitdispcoa,
                            sp_com = as.matrix(HUC07sm[,!(colnames(HUC07sm[,,i]) %in% unique(c(nnh2[["07"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC07fdis = fd_fdis(traits = HUC07traitdispcoa, sp_com = as.matrix(HUC07[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["07"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC07holdFD <- cbind(HUC07holdFD,c(HUC07fdis))
h1 <- (HUC07fdis - rowMeans(HUC07holdFD[,1:10000])) / apply(HUC07holdFD[,1:10000], 1, sd) 
HUC07$SES = h1
HUC07$fdis = HUC07fdis

saveRDS(HUC07, "./Randomizations/HUC07fdis_nng.RDS")
HUC07sm = NULL


##HUC08
HUC08sm <- readRDS("./Randomizations/HUC08sm.RDS")

HUC08traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC08[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["08"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC08traitdis <- FD::gowdis(HUC08traits)
HUC08pcoa1 <- ape::pcoa(HUC08traitdis)
HUC08traitdispcoa <- as.matrix(HUC08pcoa1$vectors)

HUC08holdFD = matrix(nrow = nrow(HUC08), ncol = 9999)
for(i in 1:9999){
  
  HUC08holdFD[,i] = fd_fdis(traits = HUC08traitdispcoa,
                            sp_com = as.matrix(HUC08sm[,!(colnames(HUC08sm[,,i]) %in% unique(c(nnh2[["08"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC08fdis = fd_fdis(traits = HUC08traitdispcoa, sp_com = as.matrix(HUC08[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["08"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC08holdFD <- cbind(HUC08holdFD,c(HUC08fdis))
h1 <- (HUC08fdis - rowMeans(HUC08holdFD[,1:10000])) / apply(HUC08holdFD[,1:10000], 1, sd) 
HUC08$SES = h1
HUC08$fdis = HUC08fdis

saveRDS(HUC08, "./Randomizations/HUC08fdis_nng.RDS")
HUC08sm = NULL

##HUC09

HUC09sm <- readRDS("./Randomizations/HUC09sm.RDS")

HUC09traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC09[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["09"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC09traitdis <- FD::gowdis(HUC09traits)
HUC09pcoa1 <- ape::pcoa(HUC09traitdis)
HUC09traitdispcoa <- as.matrix(HUC09pcoa1$vectors)

HUC09holdFD = matrix(nrow = nrow(HUC09), ncol = 9999)
for(i in 1:9999){
  
  HUC09holdFD[,i] = fd_fdis(traits = HUC09traitdispcoa,
                            sp_com = as.matrix(HUC09sm[,!(colnames(HUC09sm[,,i]) %in% unique(c(nnh2[["09"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC09fdis = fd_fdis(traits = HUC09traitdispcoa, sp_com = as.matrix(HUC09[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["09"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC09holdFD <- cbind(HUC09holdFD,c(HUC09fdis))
h1 <- (HUC09fdis - rowMeans(HUC09holdFD[,1:10000])) / apply(HUC09holdFD[,1:10000], 1, sd) 
HUC09$SES = h1
HUC09$fdis = HUC09fdis

saveRDS(HUC09, "./Randomizations/HUC09fdis_nng.RDS")
HUC09sm = NULL


##HUC10
HUC10sm <- readRDS("./Randomizations/HUC10sm.RDS")

HUC10traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC10[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["10"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC10traitdis <- FD::gowdis(HUC10traits)
HUC10pcoa1 <- ape::pcoa(HUC10traitdis)
HUC10traitdispcoa <- as.matrix(HUC10pcoa1$vectors)

HUC10holdFD = matrix(nrow = nrow(HUC10), ncol = 9999)
for(i in 1:9999){
  
  HUC10holdFD[,i] = fd_fdis(traits = HUC10traitdispcoa,
                            sp_com = as.matrix(HUC10sm[,!(colnames(HUC10sm[,,i]) %in% unique(c(nnh2[["10"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC10fdis = fd_fdis(traits = HUC10traitdispcoa, sp_com = as.matrix(HUC10[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["10"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC10holdFD <- cbind(HUC10holdFD,c(HUC10fdis))
h1 <- (HUC10fdis - rowMeans(HUC10holdFD[,1:10000])) / apply(HUC10holdFD[,1:10000], 1, sd) 
HUC10$SES = h1
HUC10$fdis = HUC10fdis

saveRDS(HUC10, "./Randomizations/HUC10fdis_nng.RDS")
HUC10sm = NULL


##HUC11
HUC11sm <- readRDS("./Randomizations/HUC11sm.RDS")

HUC11traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC11[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["11"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC11traitdis <- FD::gowdis(HUC11traits)
HUC11pcoa1 <- ape::pcoa(HUC11traitdis)
HUC11traitdispcoa <- as.matrix(HUC11pcoa1$vectors)

HUC11holdFD = matrix(nrow = nrow(HUC11), ncol = 9999)
for(i in 1:9999){
  
  HUC11holdFD[,i] = fd_fdis(traits = HUC11traitdispcoa,
                            sp_com = as.matrix(HUC11sm[,!(colnames(HUC11sm[,,i]) %in% unique(c(nnh2[["11"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC11fdis = fd_fdis(traits = HUC11traitdispcoa, sp_com = as.matrix(HUC11[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["11"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC11holdFD <- cbind(HUC11holdFD,c(HUC11fdis))
h1 <- (HUC11fdis - rowMeans(HUC11holdFD[,1:10000])) / apply(HUC11holdFD[,1:10000], 1, sd) 
HUC11$SES = h1
HUC11$fdis = HUC11fdis

saveRDS(HUC11, "./Randomizations/HUC11fdis_nng.RDS")
HUC11sm = NULL


##HUC12
HUC12sm <- readRDS("./Randomizations/HUC12sm.RDS")

HUC12traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC12[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["12"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC12traitdis <- FD::gowdis(HUC12traits)
HUC12pcoa1 <- ape::pcoa(HUC12traitdis)
HUC12traitdispcoa <- as.matrix(HUC12pcoa1$vectors)

HUC12holdFD = matrix(nrow = nrow(HUC12), ncol = 9999)
for(i in 1:9999){
  
  HUC12holdFD[,i] = fd_fdis(traits = HUC12traitdispcoa,
                            sp_com = as.matrix(HUC12sm[,!(colnames(HUC12sm[,,i]) %in% unique(c(nnh2[["12"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC12fdis = fd_fdis(traits = HUC12traitdispcoa, sp_com = as.matrix(HUC12[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["12"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC12holdFD <- cbind(HUC12holdFD,c(HUC12fdis))
h1 <- (HUC12fdis - rowMeans(HUC12holdFD[,1:10000])) / apply(HUC12holdFD[,1:10000], 1, sd) 
HUC12$SES = h1
HUC12$fdis = HUC12fdis

saveRDS(HUC12, "./Randomizations/HUC12fdis_nng.RDS")
HUC12sm = NULL


##HUC13
HUC13sm <- readRDS("./Randomizations/HUC13sm.RDS")

HUC13traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC13[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["13"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC13traitdis <- FD::gowdis(HUC13traits)
HUC13pcoa1 <- ape::pcoa(HUC13traitdis)
HUC13traitdispcoa <- as.matrix(HUC13pcoa1$vectors)

HUC13holdFD = matrix(nrow = nrow(HUC13), ncol = 9999)
for(i in 1:9999){
  
  HUC13holdFD[,i] = fd_fdis(traits = HUC13traitdispcoa,
                            sp_com = as.matrix(HUC13sm[,!(colnames(HUC13sm[,,i]) %in% unique(c(nnh2[["13"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC13fdis = fd_fdis(traits = HUC13traitdispcoa, sp_com = as.matrix(HUC13[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["13"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC13holdFD <- cbind(HUC13holdFD,c(HUC13fdis))
h1 <- (HUC13fdis - rowMeans(HUC13holdFD[,1:10000])) / apply(HUC13holdFD[,1:10000], 1, sd) 
HUC13$SES = h1
HUC13$fdis = HUC13fdis

saveRDS(HUC13, "./Randomizations/HUC13fdis_nng.RDS")
HUC13sm = NULL


##HUC14
HUC14sm <- readRDS("./Randomizations/HUC14sm.RDS")

HUC14traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC14[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["14"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC14traitdis <- FD::gowdis(HUC14traits)
HUC14pcoa1 <- ape::pcoa(HUC14traitdis)
HUC14traitdispcoa <- as.matrix(HUC14pcoa1$vectors)

HUC14holdFD = matrix(nrow = nrow(HUC14), ncol = 9999)
for(i in 1:9999){
  
  HUC14holdFD[,i] = fd_fdis(traits = HUC14traitdispcoa,
                            sp_com = as.matrix(HUC14sm[,!(colnames(HUC14sm[,,i]) %in% unique(c(nnh2[["14"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC14fdis = fd_fdis(traits = HUC14traitdispcoa, sp_com = as.matrix(HUC14[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["14"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC14holdFD <- cbind(HUC14holdFD,c(HUC14fdis))
h1 <- (HUC14fdis - rowMeans(HUC14holdFD[,1:10000])) / apply(HUC14holdFD[,1:10000], 1, sd) 
HUC14$SES = h1
HUC14$fdis = HUC14fdis

saveRDS(HUC14, "./Randomizations/HUC14fdis_nng.RDS")
HUC14sm = NULL

##HUC15
HUC15sm <- readRDS("./Randomizations/HUC15sm.RDS")

HUC15traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC15[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["15"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC15traitdis <- FD::gowdis(HUC15traits)
HUC15pcoa1 <- ape::pcoa(HUC15traitdis)
HUC15traitdispcoa <- as.matrix(HUC15pcoa1$vectors)

HUC15holdFD = matrix(nrow = nrow(HUC15), ncol = 9999)
for(i in 1:9999){
  
  HUC15holdFD[,i] = fd_fdis(traits = HUC15traitdispcoa,
                            sp_com = as.matrix(HUC15sm[,!(colnames(HUC15sm[,,i]) %in% unique(c(nnh2[["15"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC15fdis = fd_fdis(traits = HUC15traitdispcoa, sp_com = as.matrix(HUC15[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["15"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC15holdFD <- cbind(HUC15holdFD,c(HUC15fdis))
h1 <- (HUC15fdis - rowMeans(HUC15holdFD[,1:10000])) / apply(HUC15holdFD[,1:10000], 1, sd) 
HUC15$SES = h1
HUC15$fdis = HUC15fdis

saveRDS(HUC15, "./Randomizations/HUC15fdis_nng.RDS")
HUC15sm = NULL


##HUC16
HUC16sm <- readRDS("./Randomizations/HUC16sm.RDS")

HUC16traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC16[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["16"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC16traitdis <- FD::gowdis(HUC16traits)
HUC16pcoa1 <- ape::pcoa(HUC16traitdis)
HUC16traitdispcoa <- as.matrix(HUC16pcoa1$vectors)

HUC16holdFD = matrix(nrow = nrow(HUC16), ncol = 9999)
for(i in 1:9999){
  
  HUC16holdFD[,i] = fd_fdis(traits = HUC16traitdispcoa,
                            sp_com = as.matrix(HUC16sm[,!(colnames(HUC16sm[,,i]) %in% unique(c(nnh2[["16"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC16fdis = fd_fdis(traits = HUC16traitdispcoa, sp_com = as.matrix(HUC16[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["16"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC16holdFD <- cbind(HUC16holdFD,c(HUC16fdis))
h1 <- (HUC16fdis - rowMeans(HUC16holdFD[,1:10000])) / apply(HUC16holdFD[,1:10000], 1, sd) 
HUC16$SES = h1
HUC16$fdis = HUC16fdis

saveRDS(HUC16, "./Randomizations/HUC16fdis_nng.RDS")
HUC16sm = NULL


##HUC17
HUC17sm <- readRDS("./Randomizations/HUC17sm.RDS")

HUC17traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC17[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["17"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC17traitdis <- FD::gowdis(HUC17traits)
HUC17pcoa1 <- ape::pcoa(HUC17traitdis)
HUC17traitdispcoa <- as.matrix(HUC17pcoa1$vectors)

HUC17holdFD = matrix(nrow = nrow(HUC17), ncol = 9999)
for(i in 1:9999){
  
  HUC17holdFD[,i] = fd_fdis(traits = HUC17traitdispcoa,
                            sp_com = as.matrix(HUC17sm[,!(colnames(HUC17sm[,,i]) %in% unique(c(nnh2[["17"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC17fdis = fd_fdis(traits = HUC17traitdispcoa, sp_com = as.matrix(HUC17[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["17"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC17holdFD <- cbind(HUC17holdFD,c(HUC17fdis))
h1 <- (HUC17fdis - rowMeans(HUC17holdFD[,1:10000])) / apply(HUC17holdFD[,1:10000], 1, sd) 
HUC17$SES = h1
HUC17$fdis = HUC17fdis

saveRDS(HUC17, "./Randomizations/HUC17fdis_nng.RDS")
HUC17sm = NULL


##HUC18
HUC18sm <- readRDS("./Randomizations/HUC18sm.RDS")

HUC18traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC18[,-c(1:45)] %>%
                                                                dplyr::select(-any_of(c(unique(c(nnh2[["18"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC18traitdis <- FD::gowdis(HUC18traits)
HUC18pcoa1 <- ape::pcoa(HUC18traitdis)
HUC18traitdispcoa <- as.matrix(HUC18pcoa1$vectors)

HUC18holdFD = matrix(nrow = nrow(HUC18), ncol = 9999)
for(i in 1:9999){
  
  HUC18holdFD[,i] = fd_fdis(traits = HUC18traitdispcoa,
                            sp_com = as.matrix(HUC18sm[,!(colnames(HUC18sm[,,i]) %in% unique(c(nnh2[["18"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC18fdis = fd_fdis(traits = HUC18traitdispcoa, sp_com = as.matrix(HUC18[,-c(1:45)]%>%
                                                                     dplyr::select(-any_of(c(unique(c(nnh2[["18"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC18holdFD <- cbind(HUC18holdFD,c(HUC18fdis))
h1 <- (HUC18fdis - rowMeans(HUC18holdFD[,1:10000])) / apply(HUC18holdFD[,1:10000], 1, sd) 
HUC18$SES = h1
HUC18$fdis = HUC18fdis

saveRDS(HUC18, "./Randomizations/HUC18fdis_nng.RDS")
HUC18sm = NULL


##bind together

HUC01 <- readRDS("./Randomizations/HUC01fdis_nng.RDS");
HUC02 <- readRDS("./Randomizations/HUC02fdis_nng.RDS")
HUC03 <- readRDS("./Randomizations/HUC03fdis_nng.RDS");HUC04 <- readRDS("./Randomizations/HUC04fdis_nng.RDS")
HUC05 <- readRDS("./Randomizations/HUC05fdis_nng.RDS");HUC06 <- readRDS("./Randomizations/HUC06fdis_nng.RDS")
HUC07 <- readRDS("./Randomizations/HUC07fdis_nng.RDS");HUC08 <- readRDS("./Randomizations/HUC08fdis_nng.RDS")
HUC09 <- readRDS("./Randomizations/HUC09fdis_nng.RDS");HUC10 <- readRDS("./Randomizations/HUC10fdis_nng.RDS")
HUC11 <- readRDS("./Randomizations/HUC11fdis_nng.RDS");HUC12 <- readRDS("./Randomizations/HUC12fdis_nng.RDS")
HUC13 <- readRDS("./Randomizations/HUC13fdis_nng.RDS");HUC14 <- readRDS("./Randomizations/HUC14fdis_nng.RDS")
HUC15 <- readRDS("./Randomizations/HUC15fdis_nng.RDS");HUC16 <- readRDS("./Randomizations/HUC16fdis_nng.RDS")
HUC17 <- readRDS("./Randomizations/HUC17fdis_nng.RDS");HUC18 <- readRDS("./Randomizations/HUC18fdis_nng.RDS")

WHOLEHUC_nng = bind_rows(HUC01[,c(1:45,ncol(HUC01)-1,ncol(HUC01))] %>% mutate(HUC2 = "01"),
                     HUC02[,c(1:45,ncol(HUC02)-1,ncol(HUC02))]%>% mutate(HUC2 = "02"),
                     HUC03[,c(1:45,ncol(HUC03)-1,ncol(HUC03))]%>% mutate(HUC2 = "03"),
                     HUC04[,c(1:45,ncol(HUC04)-1,ncol(HUC04))]%>% mutate(HUC2 = "04"),
                     HUC05[,c(1:45,ncol(HUC05)-1,ncol(HUC05))]%>% mutate(HUC2 = "05"),
                     HUC06[,c(1:45,ncol(HUC06)-1,ncol(HUC06))]%>% mutate(HUC2 = "06"),
                     HUC07[,c(1:45,ncol(HUC07)-1,ncol(HUC07))]%>% mutate(HUC2 = "07"),
                     HUC08[,c(1:45,ncol(HUC08)-1,ncol(HUC08))]%>% mutate(HUC2 = "08"),
                     HUC09[,c(1:45,ncol(HUC09)-1,ncol(HUC09))]%>% mutate(HUC2 = "09"), 
                     HUC10[,c(1:45,ncol(HUC10)-1,ncol(HUC10))]%>% mutate(HUC2 = "10"),
                     HUC11[,c(1:45,ncol(HUC11)-1,ncol(HUC11))]%>% mutate(HUC2 = "11"),
                     HUC12[,c(1:45,ncol(HUC12)-1,ncol(HUC12))]%>% mutate(HUC2 = "12"),
                     HUC13[,c(1:45,ncol(HUC13)-1,ncol(HUC13))]%>% mutate(HUC2 = "13"),
                     HUC14[,c(1:45,ncol(HUC14)-1,ncol(HUC14))]%>% mutate(HUC2 = "14"),
                     HUC15[,c(1:45,ncol(HUC15)-1,ncol(HUC15))]%>% mutate(HUC2 = "15"), 
                     HUC16[,c(1:45,ncol(HUC16)-1,ncol(HUC16))]%>% mutate(HUC2 = "16"),
                     HUC17[,c(1:45,ncol(HUC17)-1,ncol(HUC17))]%>% mutate(HUC2 = "17"), 
                     HUC18[,c(1:45,ncol(HUC18)-1,ncol(HUC18))]%>% mutate(HUC2 = "18")
)

saveRDS(WHOLEHUC_nng,"SESNNGFinal.RDS")

##non-native####

##will need to filter species to only non-managed ones - need to do this in
##multiple spots for each HUC
gamey = read.csv("./Data/GameFishDesignation.csv")
nonnative1 = read.csv("./Data/USGSinput2_update_screened_nonnative.csv")
nonnative2 = read.csv("./Data/USGSinput1_update_screened_nonnative.csv")
nonnative = rbind(nonnative1, nonnative2)

nnh2 = nonnative %>%
  mutate(HUC2 = substr(str_pad(as.character(HUC8), 8, pad = "0"),1,2)) %>%
  dplyr::select(Scientific.Name,HUC2, Native) %>%
  filter(Native == "No") %>%
  mutate(Scientific.Name = gsub(" ","\\.",Scientific.Name)) %>%
  group_by(Scientific.Name, HUC2) %>%
  slice(1) %>%
  ungroup()

nnh2 = split(nnh2, ~HUC2)

##1,2, 3, 12-18, 4, 11, 7 are done
#10, 5, 6, 8, 9 are running

##HUC01

HUC01sm <- readRDS("./Randomizations/HUC01sm.RDS")
HUC01 = HUC01 %>%
  dplyr::select(-any_of(c("SES","fdis")))
HUC01traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC01[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["01"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC01traitdis <- FD::gowdis(HUC01traits)
HUC01pcoa1 <- ape::pcoa(HUC01traitdis)
HUC01traitdispcoa <- as.matrix(HUC01pcoa1$vectors)

HUC01holdFD = matrix(nrow = nrow(HUC01), ncol = 9999)
for(i in 1:9999){
  
  HUC01holdFD[,i] = fd_fdis(traits = HUC01traitdispcoa,
                            sp_com = as.matrix(HUC01sm[,(colnames(HUC01sm[,,i]) %in% unique(c(nnh2[["01"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC01fdis = fd_fdis(traits = HUC01traitdispcoa, sp_com = as.matrix(HUC01[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["01"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC01holdFD <- cbind(HUC01holdFD,c(HUC01fdis))
h1 <- (HUC01fdis - rowMeans(HUC01holdFD[,1:10000])) / apply(HUC01holdFD[,1:10000], 1, sd) 
HUC01$SES = h1
HUC01$fdis = HUC01fdis

saveRDS(HUC01, "./Randomizations/HUC01fdis_game.RDS")
HUC01sm = NULL

##HUC02


HUC02sm <- readRDS("./Randomizations/HUC02sm.RDS")
HUC02 = HUC02 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC02traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC02[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["02"]]$Scientific.Name, gsub(" ","\\.",gamey$Species, "Cyprinus.carpio"))),
                                                                                        "SES", "fdis")))),]

HUC02traitdis <- FD::gowdis(HUC02traits)
HUC02pcoa1 <- ape::pcoa(HUC02traitdis)
HUC02traitdispcoa <- as.matrix(HUC02pcoa1$vectors)

HUC02holdFD = matrix(nrow = nrow(HUC02), ncol = 9999)
for(i in 1:9999){
  
  HUC02holdFD[,i] = fd_fdis(traits = HUC02traitdispcoa,
                            sp_com = as.matrix(HUC02sm[,(colnames(HUC02sm[,,i]) %in% unique(c(nnh2[["02"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC02fdis = fd_fdis(traits = HUC02traitdispcoa, sp_com = as.matrix(HUC02[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["02"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC02holdFD <- cbind(HUC02holdFD,c(HUC02fdis))
h1 <- (HUC02fdis - rowMeans(HUC02holdFD[,1:10000])) / apply(HUC02holdFD[,1:10000], 1, sd) 
HUC02$SES = h1
HUC02$fdis = HUC02fdis

saveRDS(HUC02, "./Randomizations/HUC02fdis_game.RDS")
HUC02sm = NULL


##HUC03
HUC03sm <- readRDS("./Randomizations/HUC03sm.RDS")
HUC03 = HUC03 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC03traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC03[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["03"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC03traitdis <- FD::gowdis(HUC03traits)
HUC03pcoa1 <- ape::pcoa(HUC03traitdis)
HUC03traitdispcoa <- as.matrix(HUC03pcoa1$vectors)

HUC03holdFD = matrix(nrow = nrow(HUC03), ncol = 9999)
for(i in 1:9999){
  
  HUC03holdFD[,i] = fd_fdis(traits = HUC03traitdispcoa,
                            sp_com = as.matrix(HUC03sm[,(colnames(HUC03sm[,,i]) %in% unique(c(nnh2[["03"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC03fdis = fd_fdis(traits = HUC03traitdispcoa, sp_com = as.matrix(HUC03[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["03"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC03holdFD <- cbind(HUC03holdFD,c(HUC03fdis))
h1 <- (HUC03fdis - rowMeans(HUC03holdFD[,1:10000])) / apply(HUC03holdFD[,1:10000], 1, sd) 
HUC03$SES = h1
HUC03$fdis = HUC03fdis

saveRDS(HUC03, "./Randomizations/HUC03fdis_game.RDS")
HUC03sm = NULL


##HUC04
HUC04sm <- readRDS("./Randomizations/HUC04sm.RDS")
HUC04 = HUC04 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC04traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC04[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["04"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC04traitdis <- FD::gowdis(HUC04traits)
HUC04pcoa1 <- ape::pcoa(HUC04traitdis)
HUC04traitdispcoa <- as.matrix(HUC04pcoa1$vectors)

HUC04holdFD = matrix(nrow = nrow(HUC04), ncol = 9999)
for(i in 1:9999){
  
  HUC04holdFD[,i] = fd_fdis(traits = HUC04traitdispcoa,
                            sp_com = as.matrix(HUC04sm[,(colnames(HUC04sm[,,i]) %in% unique(c(nnh2[["04"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC04fdis = fd_fdis(traits = HUC04traitdispcoa, sp_com = as.matrix(HUC04[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["04"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC04holdFD <- cbind(HUC04holdFD,c(HUC04fdis))
h1 <- (HUC04fdis - rowMeans(HUC04holdFD[,1:10000])) / apply(HUC04holdFD[,1:10000], 1, sd) 
HUC04$SES = h1
HUC04$fdis = HUC04fdis

saveRDS(HUC04, "./Randomizations/HUC04fdis_game.RDS")
HUC04sm = NULL


##HUC05
HUC05sm <- readRDS("./Randomizations/HUC05sm.RDS")
HUC05 = HUC05 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC05traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC05[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["05"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC05traitdis <- FD::gowdis(HUC05traits)
HUC05pcoa1 <- ape::pcoa(HUC05traitdis)
HUC05traitdispcoa <- as.matrix(HUC05pcoa1$vectors)

HUC05holdFD = matrix(nrow = nrow(HUC05), ncol = 9999)
for(i in 1:9999){
  
  HUC05holdFD[,i] = fd_fdis(traits = HUC05traitdispcoa,
                            sp_com = as.matrix(HUC05sm[,(colnames(HUC05sm[,,i]) %in% unique(c(nnh2[["05"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC05fdis = fd_fdis(traits = HUC05traitdispcoa, sp_com = as.matrix(HUC05[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["05"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC05holdFD <- cbind(HUC05holdFD,c(HUC05fdis))
h1 <- (HUC05fdis - rowMeans(HUC05holdFD[,1:10000])) / apply(HUC05holdFD[,1:10000], 1, sd) 
HUC05$SES = h1
HUC05$fdis = HUC05fdis

saveRDS(HUC05, "./Randomizations/HUC05fdis_game.RDS")
HUC05sm = NULL

##HUC06
HUC06sm <- readRDS("./Randomizations/HUC06sm.RDS")
HUC06 = HUC06 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC06traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC06[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["06"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC06traitdis <- FD::gowdis(HUC06traits)
HUC06pcoa1 <- ape::pcoa(HUC06traitdis)
HUC06traitdispcoa <- as.matrix(HUC06pcoa1$vectors)

HUC06holdFD = matrix(nrow = nrow(HUC06), ncol = 9999)
for(i in 1:9999){
  
  HUC06holdFD[,i] = fd_fdis(traits = HUC06traitdispcoa,
                            sp_com = as.matrix(HUC06sm[,(colnames(HUC06sm[,,i]) %in% unique(c(nnh2[["06"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC06fdis = fd_fdis(traits = HUC06traitdispcoa, sp_com = as.matrix(HUC06[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["06"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC06holdFD <- cbind(HUC06holdFD,c(HUC06fdis))
h1 <- (HUC06fdis - rowMeans(HUC06holdFD[,1:10000])) / apply(HUC06holdFD[,1:10000], 1, sd) 
HUC06$SES = h1
HUC06$fdis = HUC06fdis

saveRDS(HUC06, "./Randomizations/HUC06fdis_game.RDS")
HUC06sm = NULL


##HUC07
HUC07sm <- readRDS("./Randomizations/HUC07sm.RDS")
HUC07 = HUC07 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC07traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC07[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["07"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC07traitdis <- FD::gowdis(HUC07traits)
HUC07pcoa1 <- ape::pcoa(HUC07traitdis)
HUC07traitdispcoa <- as.matrix(HUC07pcoa1$vectors)

HUC07holdFD = matrix(nrow = nrow(HUC07), ncol = 9999)
for(i in 1:9999){
  
  HUC07holdFD[,i] = fd_fdis(traits = HUC07traitdispcoa,
                            sp_com = as.matrix(HUC07sm[,(colnames(HUC07sm[,,i]) %in% unique(c(nnh2[["07"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC07fdis = fd_fdis(traits = HUC07traitdispcoa, sp_com = as.matrix(HUC07[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["07"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC07holdFD <- cbind(HUC07holdFD,c(HUC07fdis))
h1 <- (HUC07fdis - rowMeans(HUC07holdFD[,1:10000])) / apply(HUC07holdFD[,1:10000], 1, sd) 
HUC07$SES = h1
HUC07$fdis = HUC07fdis

saveRDS(HUC07, "./Randomizations/HUC07fdis_game.RDS")
HUC07sm = NULL


##HUC08
HUC08sm <- readRDS("./Randomizations/HUC08sm.RDS")
HUC08 = HUC08 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC08traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC08[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["08"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC08traitdis <- FD::gowdis(HUC08traits)
HUC08pcoa1 <- ape::pcoa(HUC08traitdis)
HUC08traitdispcoa <- as.matrix(HUC08pcoa1$vectors)

HUC08holdFD = matrix(nrow = nrow(HUC08), ncol = 9999)
for(i in 1:9999){
  
  HUC08holdFD[,i] = fd_fdis(traits = HUC08traitdispcoa,
                            sp_com = as.matrix(HUC08sm[,(colnames(HUC08sm[,,i]) %in% unique(c(nnh2[["08"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC08fdis = fd_fdis(traits = HUC08traitdispcoa, sp_com = as.matrix(HUC08[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["08"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC08holdFD <- cbind(HUC08holdFD,c(HUC08fdis))
h1 <- (HUC08fdis - rowMeans(HUC08holdFD[,1:10000])) / apply(HUC08holdFD[,1:10000], 1, sd) 
HUC08$SES = h1
HUC08$fdis = HUC08fdis

saveRDS(HUC08, "./Randomizations/HUC08fdis_game.RDS")
HUC08sm = NULL

##HUC09

HUC09sm <- readRDS("./Randomizations/HUC09sm.RDS")
HUC09 = HUC09 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC09traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC09[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["09"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC09traitdis <- FD::gowdis(HUC09traits)
HUC09pcoa1 <- ape::pcoa(HUC09traitdis)
HUC09traitdispcoa <- as.matrix(HUC09pcoa1$vectors)

HUC09holdFD = matrix(nrow = nrow(HUC09), ncol = 9999)
for(i in 1:9999){
  
  HUC09holdFD[,i] = fd_fdis(traits = HUC09traitdispcoa,
                            sp_com = as.matrix(HUC09sm[,(colnames(HUC09sm[,,i]) %in% unique(c(nnh2[["09"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC09fdis = fd_fdis(traits = HUC09traitdispcoa, sp_com = as.matrix(HUC09[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["09"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC09holdFD <- cbind(HUC09holdFD,c(HUC09fdis))
h1 <- (HUC09fdis - rowMeans(HUC09holdFD[,1:10000])) / apply(HUC09holdFD[,1:10000], 1, sd) 
HUC09$SES = h1
HUC09$fdis = HUC09fdis

saveRDS(HUC09, "./Randomizations/HUC09fdis_game.RDS")
HUC09sm = NULL


##HUC10
HUC10sm <- readRDS("./Randomizations/HUC10sm.RDS")
HUC10 = HUC10 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC10traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC10[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["10"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC10traitdis <- FD::gowdis(HUC10traits)
HUC10pcoa1 <- ape::pcoa(HUC10traitdis)
HUC10traitdispcoa <- as.matrix(HUC10pcoa1$vectors)

HUC10holdFD = matrix(nrow = nrow(HUC10), ncol = 9999)
for(i in 1:9999){
  
  HUC10holdFD[,i] = fd_fdis(traits = HUC10traitdispcoa,
                            sp_com = as.matrix(HUC10sm[,(colnames(HUC10sm[,,i]) %in% unique(c(nnh2[["10"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC10fdis = fd_fdis(traits = HUC10traitdispcoa, sp_com = as.matrix(HUC10[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["10"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC10holdFD <- cbind(HUC10holdFD,c(HUC10fdis))
h1 <- (HUC10fdis - rowMeans(HUC10holdFD[,1:10000])) / apply(HUC10holdFD[,1:10000], 1, sd) 
HUC10$SES = h1
HUC10$fdis = HUC10fdis

saveRDS(HUC10, "./Randomizations/HUC10fdis_game.RDS")
HUC10sm = NULL


##HUC11
HUC11sm <- readRDS("./Randomizations/HUC11sm.RDS")
HUC11 = HUC11 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC11traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC11[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["11"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC11traitdis <- FD::gowdis(HUC11traits)
HUC11pcoa1 <- ape::pcoa(HUC11traitdis)
HUC11traitdispcoa <- as.matrix(HUC11pcoa1$vectors)

HUC11holdFD = matrix(nrow = nrow(HUC11), ncol = 9999)
for(i in 1:9999){
  
  HUC11holdFD[,i] = fd_fdis(traits = HUC11traitdispcoa,
                            sp_com = as.matrix(HUC11sm[,(colnames(HUC11sm[,,i]) %in% unique(c(nnh2[["11"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC11fdis = fd_fdis(traits = HUC11traitdispcoa, sp_com = as.matrix(HUC11[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["11"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC11holdFD <- cbind(HUC11holdFD,c(HUC11fdis))
h1 <- (HUC11fdis - rowMeans(HUC11holdFD[,1:10000])) / apply(HUC11holdFD[,1:10000], 1, sd) 
HUC11$SES = h1
HUC11$fdis = HUC11fdis

saveRDS(HUC11, "./Randomizations/HUC11fdis_game.RDS")
HUC11sm = NULL


##HUC12
HUC12sm <- readRDS("./Randomizations/HUC12sm.RDS")
HUC12 = HUC12 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC12traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC12[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["12"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC12traitdis <- FD::gowdis(HUC12traits)
HUC12pcoa1 <- ape::pcoa(HUC12traitdis)
HUC12traitdispcoa <- as.matrix(HUC12pcoa1$vectors)

HUC12holdFD = matrix(nrow = nrow(HUC12), ncol = 9999)
for(i in 1:9999){
  
  HUC12holdFD[,i] = fd_fdis(traits = HUC12traitdispcoa,
                            sp_com = as.matrix(HUC12sm[,(colnames(HUC12sm[,,i]) %in% unique(c(nnh2[["12"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC12fdis = fd_fdis(traits = HUC12traitdispcoa, sp_com = as.matrix(HUC12[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["12"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC12holdFD <- cbind(HUC12holdFD,c(HUC12fdis))
h1 <- (HUC12fdis - rowMeans(HUC12holdFD[,1:10000])) / apply(HUC12holdFD[,1:10000], 1, sd) 
HUC12$SES = h1
HUC12$fdis = HUC12fdis

saveRDS(HUC12, "./Randomizations/HUC12fdis_game.RDS")
HUC12sm = NULL


##HUC13
HUC13sm <- readRDS("./Randomizations/HUC13sm.RDS")
HUC13 = HUC13 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC13traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC13[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["13"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC13traitdis <- FD::gowdis(HUC13traits)
HUC13pcoa1 <- ape::pcoa(HUC13traitdis)
HUC13traitdispcoa <- as.matrix(HUC13pcoa1$vectors)

HUC13holdFD = matrix(nrow = nrow(HUC13), ncol = 9999)
for(i in 1:9999){
  
  HUC13holdFD[,i] = fd_fdis(traits = HUC13traitdispcoa,
                            sp_com = as.matrix(HUC13sm[,(colnames(HUC13sm[,,i]) %in% unique(c(nnh2[["13"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC13fdis = fd_fdis(traits = HUC13traitdispcoa, sp_com = as.matrix(HUC13[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["13"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC13holdFD <- cbind(HUC13holdFD,c(HUC13fdis))
h1 <- (HUC13fdis - rowMeans(HUC13holdFD[,1:10000])) / apply(HUC13holdFD[,1:10000], 1, sd) 
HUC13$SES = h1
HUC13$fdis = HUC13fdis

saveRDS(HUC13, "./Randomizations/HUC13fdis_game.RDS")
HUC13sm = NULL


##HUC14
HUC14sm <- readRDS("./Randomizations/HUC14sm.RDS")
HUC14 = HUC14 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC14traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC14[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["14"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC14traitdis <- FD::gowdis(HUC14traits)
HUC14pcoa1 <- ape::pcoa(HUC14traitdis)
HUC14traitdispcoa <- as.matrix(HUC14pcoa1$vectors)

HUC14holdFD = matrix(nrow = nrow(HUC14), ncol = 9999)
for(i in 1:9999){
  
  HUC14holdFD[,i] = fd_fdis(traits = HUC14traitdispcoa,
                            sp_com = as.matrix(HUC14sm[,(colnames(HUC14sm[,,i]) %in% unique(c(nnh2[["14"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC14fdis = fd_fdis(traits = HUC14traitdispcoa, sp_com = as.matrix(HUC14[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["14"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC14holdFD <- cbind(HUC14holdFD,c(HUC14fdis))
h1 <- (HUC14fdis - rowMeans(HUC14holdFD[,1:10000])) / apply(HUC14holdFD[,1:10000], 1, sd) 
HUC14$SES = h1
HUC14$fdis = HUC14fdis

saveRDS(HUC14, "./Randomizations/HUC14fdis_game.RDS")
HUC14sm = NULL

##HUC15
HUC15sm <- readRDS("./Randomizations/HUC15sm.RDS")
HUC15 = HUC15 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC15traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC15[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["15"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC15traitdis <- FD::gowdis(HUC15traits)
HUC15pcoa1 <- ape::pcoa(HUC15traitdis)
HUC15traitdispcoa <- as.matrix(HUC15pcoa1$vectors)

HUC15holdFD = matrix(nrow = nrow(HUC15), ncol = 9999)
for(i in 1:9999){
  
  HUC15holdFD[,i] = fd_fdis(traits = HUC15traitdispcoa,
                            sp_com = as.matrix(HUC15sm[,(colnames(HUC15sm[,,i]) %in% unique(c(nnh2[["15"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC15fdis = fd_fdis(traits = HUC15traitdispcoa, sp_com = as.matrix(HUC15[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["15"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC15holdFD <- cbind(HUC15holdFD,c(HUC15fdis))
h1 <- (HUC15fdis - rowMeans(HUC15holdFD[,1:10000])) / apply(HUC15holdFD[,1:10000], 1, sd) 
HUC15$SES = h1
HUC15$fdis = HUC15fdis

saveRDS(HUC15, "./Randomizations/HUC15fdis_game.RDS")
HUC15sm = NULL


##HUC16
HUC16sm <- readRDS("./Randomizations/HUC16sm.RDS")
HUC16 = HUC16 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC16traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC16[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["16"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC16traitdis <- FD::gowdis(HUC16traits)
HUC16pcoa1 <- ape::pcoa(HUC16traitdis)
HUC16traitdispcoa <- as.matrix(HUC16pcoa1$vectors)

HUC16holdFD = matrix(nrow = nrow(HUC16), ncol = 9999)
for(i in 1:9999){
  
  HUC16holdFD[,i] = fd_fdis(traits = HUC16traitdispcoa,
                            sp_com = as.matrix(HUC16sm[,(colnames(HUC16sm[,,i]) %in% unique(c(nnh2[["16"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC16fdis = fd_fdis(traits = HUC16traitdispcoa, sp_com = as.matrix(HUC16[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["16"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC16holdFD <- cbind(HUC16holdFD,c(HUC16fdis))
h1 <- (HUC16fdis - rowMeans(HUC16holdFD[,1:10000])) / apply(HUC16holdFD[,1:10000], 1, sd) 
HUC16$SES = h1
HUC16$fdis = HUC16fdis

saveRDS(HUC16, "./Randomizations/HUC16fdis_game.RDS")
HUC16sm = NULL


##HUC17
HUC17sm <- readRDS("./Randomizations/HUC17sm.RDS")
HUC17 = HUC17 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC17traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC17[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["17"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC17traitdis <- FD::gowdis(HUC17traits)
HUC17pcoa1 <- ape::pcoa(HUC17traitdis)
HUC17traitdispcoa <- as.matrix(HUC17pcoa1$vectors)

HUC17holdFD = matrix(nrow = nrow(HUC17), ncol = 9999)
for(i in 1:9999){
  
  HUC17holdFD[,i] = fd_fdis(traits = HUC17traitdispcoa,
                            sp_com = as.matrix(HUC17sm[,(colnames(HUC17sm[,,i]) %in% unique(c(nnh2[["17"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC17fdis = fd_fdis(traits = HUC17traitdispcoa, sp_com = as.matrix(HUC17[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["17"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC17holdFD <- cbind(HUC17holdFD,c(HUC17fdis))
h1 <- (HUC17fdis - rowMeans(HUC17holdFD[,1:10000])) / apply(HUC17holdFD[,1:10000], 1, sd) 
HUC17$SES = h1
HUC17$fdis = HUC17fdis

saveRDS(HUC17, "./Randomizations/HUC17fdis_game.RDS")
HUC17sm = NULL


##HUC18
HUC18sm <- readRDS("./Randomizations/HUC18sm.RDS")
HUC18 = HUC18 %>%
  dplyr::select(-any_of(c("SES","fdis")))

HUC18traits <- fulldatsub[row.names(fulldatsub) %in% colnames(HUC18[,-c(1:45)] %>%
                                                                dplyr::select(any_of(c(unique(c(nnh2[["18"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                        "SES", "fdis")))),]

HUC18traitdis <- FD::gowdis(HUC18traits)
HUC18pcoa1 <- ape::pcoa(HUC18traitdis)
HUC18traitdispcoa <- as.matrix(HUC18pcoa1$vectors)

HUC18holdFD = matrix(nrow = nrow(HUC18), ncol = 9999)
for(i in 1:9999){
  
  HUC18holdFD[,i] = fd_fdis(traits = HUC18traitdispcoa,
                            sp_com = as.matrix(HUC18sm[,(colnames(HUC18sm[,,i]) %in% unique(c(nnh2[["18"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio"))),i]))$FDis
  progress(i,9999)
}

HUC18fdis = fd_fdis(traits = HUC18traitdispcoa, sp_com = as.matrix(HUC18[,-c(1:45)]%>%
                                                                     dplyr::select(any_of(c(unique(c(nnh2[["18"]]$Scientific.Name, gsub(" ","\\.",gamey$Species), "Cyprinus.carpio")),
                                                                                             "SES", "fdis")))))$FDis
HUC18holdFD <- cbind(HUC18holdFD,c(HUC18fdis))
h1 <- (HUC18fdis - rowMeans(HUC18holdFD[,1:10000])) / apply(HUC18holdFD[,1:10000], 1, sd) 
HUC18$SES = h1
HUC18$fdis = HUC18fdis

saveRDS(HUC18, "./Randomizations/HUC18fdis_game.RDS")
HUC18sm = NULL


##bind together

HUC01 <- readRDS("./Randomizations/HUC01fdis_game.RDS");
HUC02 <- readRDS("./Randomizations/HUC02fdis_game.RDS")
HUC03 <- readRDS("./Randomizations/HUC03fdis_game.RDS");HUC04 <- readRDS("./Randomizations/HUC04fdis_game.RDS")
HUC05 <- readRDS("./Randomizations/HUC05fdis_game.RDS");HUC06 <- readRDS("./Randomizations/HUC06fdis_game.RDS")
HUC07 <- readRDS("./Randomizations/HUC07fdis_game.RDS");HUC08 <- readRDS("./Randomizations/HUC08fdis_game.RDS")
HUC09 <- readRDS("./Randomizations/HUC09fdis_game.RDS");HUC10 <- readRDS("./Randomizations/HUC10fdis_game.RDS")
HUC11 <- readRDS("./Randomizations/HUC11fdis_game.RDS");HUC12 <- readRDS("./Randomizations/HUC12fdis_game.RDS")
HUC13 <- readRDS("./Randomizations/HUC13fdis_game.RDS");HUC14 <- readRDS("./Randomizations/HUC14fdis_game.RDS")
HUC15 <- readRDS("./Randomizations/HUC15fdis_game.RDS");HUC16 <- readRDS("./Randomizations/HUC16fdis_game.RDS")
HUC17 <- readRDS("./Randomizations/HUC17fdis_game.RDS");HUC18 <- readRDS("./Randomizations/HUC18fdis_game.RDS")

WHOLEHUC_game = bind_rows(HUC01[,c(1:45,ncol(HUC01)-1,ncol(HUC01))] %>% mutate(HUC2 = "01"),
                         HUC02[,c(1:45,ncol(HUC02)-1,ncol(HUC02))]%>% mutate(HUC2 = "02"),
                         HUC03[,c(1:45,ncol(HUC03)-1,ncol(HUC03))]%>% mutate(HUC2 = "03"),
                         HUC04[,c(1:45,ncol(HUC04)-1,ncol(HUC04))]%>% mutate(HUC2 = "04"),
                         HUC05[,c(1:45,ncol(HUC05)-1,ncol(HUC05))]%>% mutate(HUC2 = "05"),
                         HUC06[,c(1:45,ncol(HUC06)-1,ncol(HUC06))]%>% mutate(HUC2 = "06"),
                         HUC07[,c(1:45,ncol(HUC07)-1,ncol(HUC07))]%>% mutate(HUC2 = "07"),
                         HUC08[,c(1:45,ncol(HUC08)-1,ncol(HUC08))]%>% mutate(HUC2 = "08"),
                         HUC09[,c(1:45,ncol(HUC09)-1,ncol(HUC09))]%>% mutate(HUC2 = "09"), 
                         HUC10[,c(1:45,ncol(HUC10)-1,ncol(HUC10))]%>% mutate(HUC2 = "10"),
                         HUC11[,c(1:45,ncol(HUC11)-1,ncol(HUC11))]%>% mutate(HUC2 = "11"),
                         HUC12[,c(1:45,ncol(HUC12)-1,ncol(HUC12))]%>% mutate(HUC2 = "12"),
                         HUC13[,c(1:45,ncol(HUC13)-1,ncol(HUC13))]%>% mutate(HUC2 = "13"),
                         HUC14[,c(1:45,ncol(HUC14)-1,ncol(HUC14))]%>% mutate(HUC2 = "14"),
                         HUC15[,c(1:45,ncol(HUC15)-1,ncol(HUC15))]%>% mutate(HUC2 = "15"), 
                         HUC16[,c(1:45,ncol(HUC16)-1,ncol(HUC16))]%>% mutate(HUC2 = "16"),
                         HUC17[,c(1:45,ncol(HUC17)-1,ncol(HUC17))]%>% mutate(HUC2 = "17"), 
                         HUC18[,c(1:45,ncol(HUC18)-1,ncol(HUC18))]%>% mutate(HUC2 = "18")
)

saveRDS(WHOLEHUC_game,"SESgame_notnativeFinal.RDS")
