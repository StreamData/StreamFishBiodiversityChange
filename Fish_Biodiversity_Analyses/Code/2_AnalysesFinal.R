library(tidyverse); library(emmeans); library(ggeffects)
library(lme4); library(mgcv); library(sf);
library(vegan); library(car); library(blme);
library(adespatial); library(gratia); library(geofacet)
library(ggtern); library(glmmTMB)

emm_options(rg.limit = 1000000, lmerTest.limit = 10000, pbkrtest.limit = 10000,
            weights = "proportional")

##read in datasets from 1_WriteDataFinal
#HERE

##Stream Temperature Changes####
#summer stream temp changes
SummerstreamTemps = SummerstreamTemps %>%
  filter(wt_pred > 7)

thirds = (range(SummerstreamTemps$wt_pred)[2] - range(SummerstreamTemps$wt_pred)[1]) / 3

##means of each "stream class"; below are break points
lowerbreak = range(SummerstreamTemps$wt_pred)[1] + thirds
#15.41057
upperbreak = range(SummerstreamTemps$wt_pred)[2] - thirds
#23.80586

SummerstreamTemps = SummerstreamTemps %>%
  mutate(PastRegime = ifelse(wt_pred < lowerbreak,
                             "Cold",
                             ifelse(wt_pred > upperbreak,
                                    "Warm",
                                    "Intermediate")))

mc = lm(estimate*10 ~ PastRegime, data = SummerstreamTemps)
Anova(mc,test = "F")
pairs(emmeans(mc, ~PastRegime), adjust = "none")
# contrast            estimate      SE   df t.ratio p.value
# Cold - Intermediate   0.1085 0.00717 3399  15.138  <.0001
# Cold - Warm           0.1788 0.00788 3399  22.690  <.0001
# Intermediate - Warm   0.0703 0.00490 3399  14.331  <.0001
test(emmeans(mc, ~PastRegime), adjust = "none")
# PastRegime   emmean      SE   df t.ratio p.value
# Cold         0.2619 0.00669 3399  39.152  <.0001
# Intermediate 0.1533 0.00258 3399  59.341  <.0001
# Warm         0.0831 0.00417 3399  19.937  <.0001

##Diversity analyses####
###Whole community####

#CPUE
m1 <- lmer(log(TotalCPUE100m + 0.01) ~ 
             Year*HUC2 + Agency +
             StreamOrder  +
             SANDCAT+
             SampleTypeCode + log(MethodEffort)+ log(ReachLengthFished_100m)+
             log(PredictedWettedWidth_m) +
             poly(log(WholeConductivity), 2, raw = TRUE) +
             Year*wt_pred_new +
             Landuse*Year+
             (1|SiteNumber) + (1|YearC),
           data = fishdatCPUE)

Anova(m1, type = 3, test = "F")

performance::r2(m1)
#  Conditional R2: 0.741
#  Marginal R2: 0.450

test(emtrends(m1, ~wt_pred_new, var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))

emtrends(m1, ~Landuse,
         weights = "proportional", var = "Year")

pairs(emtrends(m1, ~Landuse,
               weights = "proportional", var = "Year"), adjust = "none")
# contrast                       estimate      SE   df t.ratio p.value
# Agriculture - (Forest/Wetland) -0.02138 0.00882 4141  -2.423  0.0154
# Agriculture - Grassland         0.00294 0.00967 4163   0.304  0.7612
# Agriculture - Urban            -0.00835 0.00909 4142  -0.919  0.3582
# (Forest/Wetland) - Grassland    0.02432 0.00807 4079   3.012  0.0026
# (Forest/Wetland) - Urban        0.01303 0.00725 3976   1.796  0.0726
# Grassland - Urban              -0.01129 0.00876 3709  -1.288  0.1978

#rarefied richness
m_rare <- lmer(log(RarefiedRichness) ~ Year*HUC2+Agency +
                 StreamOrder + SampleTypeCode +
                 SANDCAT+
                 log(PredictedWettedWidth_m) +
                 poly(log(WholeConductivity),2,raw = TRUE) +
                 Year*wt_pred_new +
                 Landuse+
                 (1|SiteNumber) + (1|YearC),
               data = fishdat_rarefiedRich,
               control = lmerControl(optimizer="nloptwrap",optCtrl=list(maxfun=2e7,
                                                                     xtol_abs=1e-8,
                                                                     ftol_abs=1e-8)))



car::Anova(m_rare, type = 3, test = "F")

performance::r2(m_rare)
# Conditional R2: 0.879
# Marginal R2: 0.534

test(emtrends(m_rare, ~wt_pred_new, var = "Year", 
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))

emmeans(m_rare, ~ Landuse,
        weights = "proportional")

pairs(emmeans(m_rare, ~ Landuse,
              weights = "proportional"), adjust = "none")
# contrast                       estimate     SE   df t.ratio p.value
# Agriculture - (Forest/Wetland)  0.00536 0.0268 2679   0.200  0.8414
# Agriculture - Grassland         0.03940 0.0286 2717   1.379  0.1680
# Agriculture - Urban             0.12642 0.0335 2515   3.772  0.0002
# (Forest/Wetland) - Grassland    0.03404 0.0245 2777   1.387  0.1656
# (Forest/Wetland) - Urban        0.12105 0.0285 2558   4.240  <.0001
# Grassland - Urban               0.08701 0.0328 2537   2.653  0.0080

#SES
m1f <- blmer(SES ~ Year*HUC2 + Agency + StreamOrder  + SampleTypeCode +
               SANDCAT+
               log(PredictedWettedWidth_km)+
               poly(log(WholeConductivity), 2, raw = TRUE)+
               wt_pred_new*Year+
               Landuse+
               (1|SiteNumber) + (1|YearC),
             data = fishdat,
             control = lmerControl(optimizer="Nelder_Mead"))

performance::r2(m1f)
# Conditional R2: 0.541
# Marginal R2: 0.126

Anova(m1f, type = 3, test = "F")

test(emtrends(m1f,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))

emmeans(m1f, ~ Landuse,
        weights = "proportional")

pairs(emmeans(m1f, ~ Landuse,
              weights = "proportional"), adjust = "none")
# contrast                       estimate     SE   df t.ratio p.value
# Agriculture - (Forest/Wetland)   0.1838 0.0624 2787   2.946  0.0032
# Agriculture - Grassland          0.1506 0.0657 2722   2.293  0.0219
# Agriculture - Urban              0.0848 0.0767 2534   1.105  0.2692
# (Forest/Wetland) - Grassland    -0.0332 0.0582 2799  -0.571  0.5682
# (Forest/Wetland) - Urban        -0.0990 0.0667 2584  -1.484  0.1378
# Grassland - Urban               -0.0657 0.0754 2499  -0.872  0.3835

#LCBD
mLCBD <- lmer(LCBD_s2 ~ Year*HUC2+Agency + StreamOrder  + SampleTypeCode +
                log(PredictedWettedWidth_km)+
                SANDCAT+
                poly(log(WholeConductivity), 2, raw = TRUE)+
                Year * wt_pred_new +
                Landuse+
                (1|SiteNumber) + (1|YearC),
              data = fishdat_LCBD)

Anova(mLCBD, type = 3, test = "F")

performance::r2(mLCBD)
# Conditional R2: 0.777
# Marginal R2: 0.087

test(emtrends(mLCBD,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))

emmeans(mLCBD, ~ Landuse,
        weights = "proportional")

pairs(emmeans(mLCBD, ~ Landuse,
              weights = "proportional"), adjust = "none")
# contrast                       estimate      SE   df t.ratio p.value
# Agriculture - (Forest/Wetland) -0.02501 0.00778 3191  -3.213  0.0013
# Agriculture - Grassland        -0.01886 0.00815 3189  -2.314  0.0207
# Agriculture - Urban            -0.01369 0.00977 3008  -1.402  0.1611
# (Forest/Wetland) - Grassland    0.00615 0.00695 3303   0.884  0.3765
# (Forest/Wetland) - Urban        0.01131 0.00842 3037   1.343  0.1794
# Grassland - Urban               0.00516 0.00944 3058   0.547  0.5846

###native, non-game####
#CPUE
m1_native <- lmer(log(TotalCPUE100m + 0.01) ~ Year*HUC2+Agency + StreamOrder  +
                    SampleTypeCode + log(MethodEffort)+ log(ReachLengthFished_100m)+
                    SANDCAT+log(PredictedWettedWidth_m) +
                    poly(log(WholeConductivity), 2, raw = TRUE) +
                    Landuse*Year +
                    Year*wt_pred_new + 
                    (1|SiteNumber) + (1|YearC),
                  data = fishdatCPUE_native)

Anova(m1_native, type = 3, test = "F")

performance::r2(m1_native)

test(emtrends(m1_native,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))

#richness
m_rare_native<- lmer(log(RarefiedRichness) ~ Year*HUC2+Agency +
                       StreamOrder + SampleTypeCode +
                       SANDCAT+
                       log(PredictedWettedWidth_m) +
                       poly(log(WholeConductivity),2,raw = TRUE) +
                       Year*wt_pred_new +
                       Landuse+
                       (1|SiteNumber) + (1|YearC),
                     data = fishdat_rarefiedRich_nng,
                     control = lmerControl(optimizer = "nlminbwrap"))

car::Anova(m_rare_native, type = 3, test = "F")

performance::r2(m_rare_native)

test(emtrends(m_rare_native, ~wt_pred_new, var = "Year", pbkrtest.limit = 60000,
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))

#SES
m1f_native <- blmer(SES ~ Year*HUC2 + Agency + StreamOrder  + SampleTypeCode +
                   SANDCAT+
                   log(PredictedWettedWidth_km)+
                   poly(log(WholeConductivity), 2, raw = TRUE)+
                   wt_pred_new*Year+
                   Landuse+
                   (1|SiteNumber) + (1|YearC),
                 data = fishdat_native,
                 control = lmerControl(optimizer = "nlminbwrap"))

summary(m1f_native)
performance::r2(m1f_native)
Anova(m1f_native, type = 3, test = "F")

test(emtrends(m1f_native,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))

#LCBD
mLCBD_native <- lmer(LCBD2_s2 ~ Year*HUC2+Agency + StreamOrder  + SampleTypeCode +
                       log(PredictedWettedWidth_km)+
                       SANDCAT+
                       poly(log(WholeConductivity), 2, raw = TRUE)+
                       Year * wt_pred_new +
                       Landuse+
                       (1|SiteNumber) + (1|YearC),
                     data = fishdat_native_LCBD,
                     control = lmerControl(optimizer = "nlminbwrap"))

performance::r2(mLCBD_native)

Anova(mLCBD_native, type = 3, test = "F")

test(emtrends(mLCBD_native,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))



###game, non-native####
#CPUE
#CPUE
m1_notnative <- lmer(log(TotalCPUE100m + 0.01) ~ Year*HUC2+Agency + StreamOrder  +
                    SampleTypeCode + log(MethodEffort)+ log(ReachLengthFished_100m)+
                    SANDCAT+log(PredictedWettedWidth_m) +
                    poly(log(WholeConductivity), 2, raw = TRUE) +
                    Landuse +
                    Year*wt_pred_new + 
                    (1|SiteNumber) + (1|YearC),
                  data = fishdatCPUE_notnative %>% filter(TotalCPUE100m != 0))

Anova(m1_notnative, type = 3, test = "F")

performance::r2(m1_notnative)

test(emtrends(m1_notnative,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))


#richness
m_rare_notnative<- lmer(log(RarefiedRichness) ~ Year*HUC2+Agency +
                       StreamOrder + SampleTypeCode +
                       SANDCAT+
                       log(PredictedWettedWidth_m) +
                       poly(log(WholeConductivity),2,raw = TRUE) +
                       Year*wt_pred_new +
                       Landuse+
                       (1|SiteNumber) + (1|YearC),
                     data = fishdat_rarefiedRich_notnative)

Anova(m_rare_notnative, type = 3, test = "F")

performance::r2(m_rare_notnative)

test(emtrends(m_rare_notnative, ~wt_pred_new, var = "Year", pbkrtest.limit = 60000,
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))

##ses
m1f_notnative<- blmer(SES ~ Year*HUC2 + Agency + StreamOrder  + SampleTypeCode +
                   SANDCAT+
                   log(PredictedWettedWidth_km)+
                   poly(log(WholeConductivity), 2, raw = TRUE)+
                   wt_pred_new*Year+
                   Landuse+
                   (1|SiteNumber) + (1|YearC),
                 data = fishdat_notnativeSES,
                 control = lmerControl(optimizer = "nlminbwrap"))

summary(m1f_notnative)
performance::r2(m1f_notnative)
Anova(m1f_notnative,type = 3, test = "F")

test(emtrends(m1f_notnative,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"
              ))

##lcbd
mLCBD_notnative <- lmer(LCBD2_s2 ~ Year*HUC2+Agency + StreamOrder  + SampleTypeCode +
                       log(PredictedWettedWidth_km)+
                       SANDCAT+
                       poly(log(WholeConductivity), 2, raw = TRUE)+
                       Year * wt_pred_new +
                       Landuse+
                       (1|SiteNumber) + (1|YearC),
                     data = fishdat_notnative_LCBD)

Anova(mLCBD_notnative, type = 3, test = "F")

performance::r2(mLCBD_notnative)

test(emtrends(mLCBD_notnative, ~wt_pred_new, var = "Year", pbkrtest.limit = 60000,
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))

##OPE####
###trends####
mOPE <- glmmTMB(Opportunist ~ Year*HUC2+
                  Agency + 
                  Year * wt_pred_new +
                  Landuse+
                  (1|SiteNumber) +
                  (1|YearC),
                family = beta_family(link = "logit"),
                data = fishdatOPE,
                control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

summary(mOPE)
Anova(mOPE, type = 3)

test(emtrends(mOPE, ~ wt_pred_new, var = "Year",
              weights = "proportional",
              at = list(wt_pred_new = c(13.0,20.5,25.6))))


mOPE2 <- glmmTMB(Periodic ~ Year*HUC2+
                   Agency + 
                   Year * wt_pred_new +
                   Landuse+
                   (1|SiteNumber) +
                   (1|YearC),
                 family = beta_family(link = "logit"),
                 data = fishdatOPE,
                 control=glmmTMBControl(optimizer=nlminb,
                                        optCtrl=list(iter.max=1e3,
                                                     eval.max=1e3)))

summary(mOPE2)
Anova(mOPE2, type = 3)

test(emtrends(mOPE2, ~ wt_pred_new, var = "Year",
              weights = "proportional",
              at = list(wt_pred_new = c(13.0,20.5,25.6))))

mOPE3 <- glmmTMB(Equilibrium ~ Year*HUC2+
                   Agency + 
                   Year * wt_pred_new +
                   Landuse+
                   (1|SiteNumber) +
                   (1|YearC),
                 family = beta_family(link = "logit"),
                 data = fishdatOPE,
                 control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

Anova(mOPE3, type = 3)





##change in local fish related to Climate change and introduced/game----

##test for local trends predicted by change in temperature (from lm() above),
##change in managed (introduced or game)
##
##gather site-level estimates by getting the year coefficient (temporal trend),
##for each unique combination of site-level covariates that interact with year
##in the diversity analyses (e.g., HUC, landuse, and past stream temperature)
##and propagate the error
##Do this for all local and non-game diversity endpoints as well as  
##introduced or game abundance

emm_options(rg.limit = 100000000, lmerTest.limit = 10000, pbkrtest.limit = 10000,
            weights = "proportional")
# 
# CPUEtrendswhole = emtrends(m1_native, ~wt_pred_new|HUC2|Landuse, var = "Year",
#                            at = list(wt_pred_new = unique(fishdatCPUE$wt_pred_new)),
#                            lmer.df = "asymptotic",
#                            weights = "proportional")
# saveRDS(CPUEtrendswhole, "./Data/CPUETrends.RDS")
CPUEtrendswhole = readRDS("./Data/CPUETrends.RDS")

CPUEtrendswhole = data.frame(CPUEtrendswhole) %>%
  rename(CPUE.trend = Year.trend,
         CPUE.trend.SE = SE)


# Richtrendswhole = emtrends(m_rare_native, ~wt_pred_new|HUC2, var = "Year",
#                            at = list(wt_pred_new = unique(fishdatCPUE$wt_pred_new)),
#                            lmer.df = "asymptotic",
#                            weights = "proportional")
# saveRDS(Richtrendswhole, "./Data/RichTrends.RDS")
Richtrendswhole = readRDS("./Data/RichTrends.RDS")

Richtrendswhole = data.frame(Richtrendswhole) %>%
  rename(Rich.trend = Year.trend,
         Rich.trend.SE = SE)

# LCBDtrendswhole = emtrends(mLCBD_native, ~wt_pred_new|HUC2, var = "Year",
#                            at = list(wt_pred_new = unique(fishdat_LCBD$wt_pred_new)),
#                            lmer.df = "asymptotic",
#                            weights = "proportional")
# saveRDS(LCBDtrendswhole, "./Data/LCBDTrends.RDS")
LCBDtrendswhole = readRDS("./Data/LCBDTrends.RDS")

LCBDtrendswhole = data.frame(LCBDtrendswhole) %>%
  rename(LCBD.trend = Year.trend,
         LCBD.trend.SE = SE)

# SEStrendswhole = emtrends(m1f_native, ~wt_pred_new|HUC2, var = "Year",
#                           at = list(wt_pred_new = unique(fishdat_LCBD$wt_pred_new)),
#                           lmer.df = "asymptotic",
#                           weights = "proportional")
# saveRDS(SEStrendswhole, "./Data/SESTrends.RDS")
SEStrendswhole = readRDS("./Data/SESTrends.RDS")

SEStrendswhole = data.frame(SEStrendswhole) %>%
  rename(SES.trend = Year.trend,
         SES.trend.SE = SE)




# CPUEtrendsintgame = emtrends(m1_notnative, ~wt_pred_new|HUC2, var = "Year",
#                              at = list(wt_pred_new = unique(fishdatCPUE$wt_pred_new)),
#                              lmer.df = "asymptotic",
#                              weights = "proportional")
# saveRDS(CPUEtrendsintgame, "./Data/CPUETrendsintgame.RDS")
CPUEtrendsintgame = readRDS("./Data/CPUETrendsintgame.RDS")

CPUEtrendsintgame = data.frame(CPUEtrendsintgame) %>%
  rename(intgame.CPUE.trend = Year.trend,
         intgame.CPUE.trend.SE = SE)

trendies = fishdatCPUE %>%
  group_by(COMID) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(COMID, wt_pred_new, Landuse, HUC2) %>%
  left_join(CPUEtrendswhole %>% dplyr::select(-df, -asymp.LCL, -asymp.UCL))%>%
  left_join(Richtrendswhole %>% dplyr::select(-df, -asymp.LCL, -asymp.UCL)) %>%
  left_join(LCBDtrendswhole %>% dplyr::select(-df, -asymp.LCL, -asymp.UCL)) %>%
  left_join(SEStrendswhole %>% dplyr::select(-df, -asymp.LCL, -asymp.UCL)) %>%
  left_join(CPUEtrendsintgame %>% dplyr::select(-df, -asymp.LCL, -asymp.UCL)) %>%
  left_join(SummerstreamTemps %>% dplyr::select(-wt_pred, -p.value)) %>%
  filter(estimate > -0.05)

trendies$CPUE_vi = trendies$CPUE.trend.SE^2
trendies$Rich_vi = trendies$Rich.trend.SE^2
trendies$LCBD_vi = trendies$LCBD.trend.SE^2
trendies$SES_vi = trendies$SES.trend.SE^2

##scaled version first to get standardized coefficients
##then normal to generate figures
CPUEtrends = lm(scale(CPUE.trend) ~ scale(estimate)*scale(intgame.CPUE.trend),
                data = trendies,
                weights = 1/CPUE_vi)
summary(CPUEtrends)

CPUEtrends = lm(CPUE.trend ~ estimate*intgame.CPUE.trend,
                data = trendies,
                weights = 1/CPUE_vi)
Anova(CPUEtrends, type = 3)

##
Richtrends = lm(scale(Rich.trend) ~ scale(estimate)*scale(intgame.CPUE.trend),
                data = trendies,
                weights = 1/Rich_vi)
summary(Richtrends)

Richtrends = lm(Rich.trend ~ estimate*intgame.CPUE.trend,
                data = trendies,
                weights = 1/Rich_vi)
Anova(Richtrends, type = 3)


##
LCBDtrends = lm(scale(LCBD.trend) ~ scale(estimate)*scale(intgame.CPUE.trend),
                data = trendies,
                weights = 1/LCBD_vi)
summary(LCBDtrends)

LCBDtrends = lm(LCBD.trend ~ estimate*intgame.CPUE.trend,
                data = trendies,
                weights = 1/LCBD_vi)
Anova(LCBDtrends, type = 3)

##
SEStrends = lm(scale(SES.trend) ~ scale(estimate)*scale(intgame.CPUE.trend),
               data = trendies,
               weights = 1/SES_vi)
summary(SEStrends)

SEStrends = lm(SES.trend ~ estimate*intgame.CPUE.trend,
               data = trendies,
               weights = 1/SES_vi)
Anova(SEStrends, type = 3)


