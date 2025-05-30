library(tidyverse);library(emmeans);library(ggeffects);
library(cowplot); library(scales)
library(tmap); library(sf); library(nngeo)
#Main text figures----
##figure 2 ----
SummerstreamTemps = SummerstreamTemps %>%
  filter(wt_pred > 7)

thirds = (range(SummerstreamTemps$wt_pred)[2] - range(SummerstreamTemps$wt_pred)[1]) / 3

##means of each "stream class"; below are break points
lowerbreak = range(SummerstreamTemps$wt_pred)[1] + thirds
#15.41057
upperbreak = range(SummerstreamTemps$wt_pred)[2] - thirds
#23.80586

JulystreamTemps2 = SummerstreamTemps %>%
  group_by(COMID) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(PastRegime = ifelse(wt_pred < lowerbreak,
                              "Cold",
                              ifelse(wt_pred > upperbreak,
                                     "Warm",
                                     "Intermediate")))
##mean values
JulystreamTemps2 %>%
  group_by(PastRegime) %>%
  summarize(mean = mean(wt_pred),
            count =n()) %>%
  ungroup() %>%
  mutate(prop = count / sum(count))


hist1 = ggplot(JulystreamTemps2, aes(x = wt_pred,
                                     fill = PastRegime))+
  # facet_wrap(~CollectionYear, scales = "free_y")+
  geom_histogram(binwidth = 0.5)+
  geom_vline(xintercept = 13.0, linetype = "dashed")+
  geom_vline(xintercept = 20.5, linetype = "dashed")+
  geom_vline(xintercept = 25.6, linetype = "dashed")+
  annotate("text", label = "13.0",
           y = 215, x = 10.9, size = 2.5)+
  annotate("text", label = "20.5",
           y = 215, x = 18.25, size = 2.5)+
  annotate("text", label = "25.6",
           y = 215, x = 28, size = 2.5)+
  # xlab(expression("Mean summer temperature " (degree*C*", 1990-94")))+
  xlab("Mean summer temperature (1990-94)")+
  ylab("Number of sites")+
  scale_y_continuous(limits = c(0,225),
                     expand = c(0,0))+
  scale_fill_manual(name = "Past temperature regime",
                    values = c("blue","violet","red"))+
  scale_x_continuous(labels = scales::label_number(suffix = "°C"))+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.margin = margin(c(-5,0,0,0)),
        axis.text = element_text(color = "black",
                                   size = 8),
        axis.title = element_text(size = 9),
        axis.title.x = element_text(hjust = 1),
        plot.margin = margin(r = 5))

# legendA = get_plot_component(hist1 + theme(legend.box.margin = margin(c(0,0,0,0))), 'guide-box-bottom', return_all = TRUE)
hist12 <- hist1 + theme(legend.position='none')

#MAP


HUC2 <- read_sf("./Data/HUC2_ClipUS.shp")
HUC2 <- st_transform(HUC2,crs = 5070)

US_Bounds <- read_sf("./Data/UnitedStates.shp")
US_Bounds <- st_transform(US_Bounds,crs = 5070)

pts <- readRDS("./Data/sitesLatLong 3.rds") %>%
  mutate(COMID = as.character(COMID))
pts = pts %>% left_join(SummerstreamTemps, by = join_by(COMID == COMID))
pts <- pts[!is.na(pts$wt_pred),]

pts <- st_as_sf(pts, coords = c("Longitude_dd","Latitude_dd"), crs = 4269)
pts <- st_transform(pts,crs = 5070)
pts$Agency <- factor(pts$Agency)

#< 15.817 is cold; >23.215 is warm
pts$Temp_Class <- ifelse(pts$wt_pred <= 18, "Cold", 
                         ifelse(pts$wt_pred >= 22, "Warm", 
                                "Intermediate"))
pts$Temp_Class <- factor(pts$Temp_Class, 
                         levels = c("Cold", "Intermediate","Warm"))

Plot_Thermal <- tm_shape(HUC2)+
  tm_borders()+
  tm_shape(US_Bounds)+
  tm_borders()+
  tm_shape(pts)+
  tm_dots(col = "Temp_Class", 
          palette = c("blue", "violet", "red"),
          size = 0.125, shape = 21, legend.show = F) +
  # tm_add_legend(type = "symbol", size=0.5,
  #               shape = 21, col = c("blue", "violet", "red"),
  #               labels = c("Cold", "Intermediate","Warm"),
  #               title = "Past temperature regime", 
  #               is.portrait = T)+
  tm_scale_bar(position = c(0,0),
               breaks = c(0,500,1000),
               text.size = .5, lwd = 1, width = .25)+
  tm_layout(frame = FALSE, 
            inner.margins = c(0,0.0,0.0,0),
            asp = 1.6,
            bg.color = "transparent")

##make sure plot window is sized so that ratio is 1.1
MAPHOLD = tmap_grob(Plot_Thermal)


#add legend to bottom center of the figure
##fake code just for a legend


hld1 = ggplot(JulystreamTemps2, aes(x = wt_pred,
                                    y = estimate,
                                     fill = PastRegime))+
  # facet_wrap(~CollectionYear, scales = "free_y")+
geom_point(shape = 21, color = "grey43")+
  scale_fill_manual(name = "Past temperature regime",
                    values = c("blue","violet","red"))+
  theme_bw()+
  theme(legend.position = "bottom")

legend.A = (get_plot_component(hld1 + theme(legend.margin = margin(-10,0,-10,0)), 'guide-box-bottom', return_all = TRUE))


fullfig1 = plot_grid(plot_grid(NULL,hist12,NULL,
                               nrow = 3,
                               rel_heights = c(0.25,1,0.25),
                               labels = c("","A",""),
                               label_y = 1.1),
                     NULL,
                     MAPHOLD,
                     NULL,
          nrow = 1,
          labels = c("","","B",""),
          rel_widths = c(.8,.05,1,.05),
          label_x = -.05,
          label_y = .9)

fullfig1.1 = plot_grid(fullfig1,NULL,legend.A,
                       ncol=1,
                       nrow = 3,
                       rel_heights = c(1,-.2,.2))

ggsave("./Figures/Figure1AB_new.jpg",fullfig1.1,
       height = 3.25, width = 5.25, units = "in",
       dpi = 600)
ggsave("./Figures/Figure1AB_new.eps",fullfig1.1,
       height = 3.25, width = 5.25, units = "in",
       dpi = 600)

##Figure 3 Past stream temp forest plot######
#CPUE
test(emtrends(m1, ~wt_pred_new, var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE    df t.ratio p.value
# 13.0  -0.028292 0.00912 563.8  -3.102  0.0020
# 20.5   0.000305 0.00462  44.8   0.066  0.9477
# 25.6   0.019751 0.00686 227.8   2.879  0.0044

test(emtrends(m1_native,  ~ wt_pred_new,
              var = "Year", at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE    df t.ratio p.value
# 13.0   -0.05935 0.01222 673.1  -4.856  <.0001
# 20.5   -0.00345 0.00588  47.7  -0.587  0.5600
# 25.6    0.03456 0.00914 298.2   3.783  0.0002

test(emtrends(m1.sp,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE  df z.ratio p.value
# 13.0   -0.02610 0.00887 Inf  -2.942  0.0033
# 20.5   -0.00161 0.00433 Inf  -0.373  0.7092
# 25.6    0.01504 0.00663 Inf   2.267  0.0234

test(emtrends(m1_notnative,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE    df t.ratio p.value
# 13.0    0.00350 0.01146 500.2   0.305  0.7605
# 20.5    0.00871 0.00591  44.7   1.475  0.1473
# 25.6    0.01226 0.00841 199.9   1.457  0.1468

CPUEtrends = data.frame(Year.trend = c(-0.028292,0.000305,0.019751,
                                       -0.05935,-0.00345,0.03456,
                                       -0.02610,-0.00161,0.01504,
                                       0.00350,0.00871,0.01226),
                        SE = c(0.00912,0.00462,0.00686,
                               0.01222,0.00588,0.00914,
                               0.00887,0.00433,0.00663,
                               0.01146,0.00591,0.00841),
                        Temperature = rep(c("Cold","Intermediate","Warm"), times = 4),
                        Type = rep(c("Whole community","Native, non-game","Spatial model","Non-native, game"), each = 3),
                        Endpoint = "Abundance")

CPUEtrends$decadepercchange = exp(CPUEtrends$Year.trend)^10 - 1


#richness
test(emtrends(m_rare, ~wt_pred_new, var = "Year", pbkrtest.limit = 60000,
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE    df t.ratio p.value
# 13.0   -0.01428 0.00288 866.6  -4.962  <.0001
# 20.5   -0.00258 0.00123  42.7  -2.100  0.0417
# 25.6    0.00538 0.00197 316.1   2.729  0.0067

test(emtrends(m_rare_native, ~wt_pred_new, var = "Year", pbkrtest.limit = 60000,
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE     df t.ratio p.value
# 13.0   -0.00920 0.00382 1483.7  -2.408  0.0162
# 20.5   -0.00278 0.00125   33.3  -2.225  0.0329
# 25.6    0.00159 0.00221  349.6   0.720  0.4723


test(emtrends(m_rare.sp,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE  df z.ratio p.value
# 13.0   -0.01367 0.00282 Inf  -4.851  <.0001
# 20.5   -0.00270 0.00110 Inf  -2.457  0.0140
# 25.6    0.00476 0.00190 Inf   2.499  0.0125

test(emtrends(m_rare_notnative, ~wt_pred_new, var = "Year", pbkrtest.limit = 60000,
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE    df t.ratio p.value
# 13.0   -0.01368 0.00377 637.2  -3.625  0.0003
# 20.5   -0.00396 0.00183  46.3  -2.163  0.0357
# 25.6    0.00264 0.00258 182.8   1.024  0.3071

Richtrends = data.frame(Year.trend = c(-0.01428,-0.00258,0.00538,
                                       -0.00920,-0.00278,0.00159,
                                       -0.01367,-0.00270,0.00476,
                                       -0.01368,-0.00396,0.00264),
                        SE = c(0.00288,0.00123,0.00197,
                               0.00382,0.00125,0.00221,
                               0.00282,0.00110,0.00190,
                               0.00377,0.00183,0.00258),
                        Temperature = rep(c("Cold","Intermediate","Warm"), times = 4),
                        Type = rep(c("Whole community","Native, non-game","Spatial model","Non-native, game"), each = 3),
                        Endpoint = "Rarefied species richness")
Richtrends$decadepercchange = exp(Richtrends$Year.trend)^10 - 1

#SES

test(emtrends(m1f,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE    df t.ratio p.value
# 13.0    0.01009 0.00802 863.9   1.258  0.2088
# 20.5    0.00337 0.00352  41.7   0.956  0.3445
# 25.6   -0.00120 0.00530 234.7  -0.227  0.8206

test(emtrends(m1f_native,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE  df t.ratio p.value
# 13.0   0.000707 0.00750 984   0.094  0.9250
# 20.5  -0.001720 0.00301  36  -0.572  0.5712
# 25.6  -0.003370 0.00486 296  -0.693  0.4887

test(emtrends(m1f.sp,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE  df z.ratio p.value
# 13.0    0.01556 0.00739 Inf   2.105  0.0353
# 20.5    0.00266 0.00271 Inf   0.983  0.3255
# 25.6   -0.00611 0.00491 Inf  -1.245  0.2132

test(emtrends(m1f_notnative, ~wt_pred_new,
              var = "Year", 
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE    df t.ratio p.value
# 13.0    0.01404 0.00817 923.5   1.719  0.0860
# 20.5   -0.00125 0.00340  36.8  -0.369  0.7145
# 25.6   -0.01166 0.00532 256.2  -2.192  0.0293

fdistrends = data.frame(Year.trend = c(0.01009,0.00337,-0.00120,
                                       -0.000707,-0.001720,-0.003370,
                                       0.01556,0.00266,-0.00611,
                                       0.01404,-0.00125, -0.01166),
                        SE = c(0.00802,0.00352,0.00530,
                               0.00750,0.00301,0.00486,
                               0.00739,0.00271,0.00491,
                               0.00817,0.00340,0.00532),
                        Temperature = rep(c("Cold","Intermediate","Warm"), times = 4),
                        Type = rep(c("Whole community","Native, non-game","Spatial model","Non-native, game"), each = 3),
                        Endpoint = "Functional diversity")

#LCBD

test(emtrends(mLCBD,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend       SE    df t.ratio p.value
# 13.0   0.003123 0.000814 542.9   3.835  0.0001
# 20.5  -0.000506 0.000407  45.3  -1.245  0.2195
# 25.6  -0.002974 0.000585 212.7  -5.080  <.0001

test(emtrends(mLCBD_native,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend       SE     df t.ratio p.value
# 13.0   0.000537 0.001151 1001.6   0.466  0.6410
# 20.5  -0.000463 0.000491   49.7  -0.944  0.3499
# 25.6  -0.001144 0.000718  244.5  -1.594  0.1123

test(emtrends(mLCBD.sp,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE  df z.ratio p.value
# 13.0   0.002672 0.00149 Inf   1.791  0.0732
# 20.5  -0.000294 0.00132 Inf  -0.222  0.8245
# 25.6  -0.002310 0.00139 Inf  -1.668  0.0954

test(emtrends(mLCBD_notnative, ~wt_pred_new, var = "Year", pbkrtest.limit = 60000,
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend       SE    df t.ratio p.value
# 13.0    0.00262 0.000945 652.0   2.773  0.0057
# 20.5   -0.00113 0.000401  37.5  -2.814  0.0077
# 25.6   -0.00368 0.000648 357.3  -5.683  <.0001

LCBDtrends = data.frame(Year.trend = c(0.003123,-0.000506,-0.002974,
                                       0.000537,-0.000463,-0.001144,
                                       0.002672,-0.000294,-0.002310,
                                       0.00262,-0.00113,-0.00368),
                        SE = c(0.000814,0.000407,0.000585,
                               0.001151,0.000491,0.000718,
                               0.00149, 0.00132, 0.00139,
                               0.000945,0.000401,0.000648
                               ),
                        Temperature = rep(c("Cold","Intermediate","Warm"), times = 4),
                        Type = rep(c("Whole community","Native, non-game","Spatial model","Non-native, game"), each = 3),
                        Endpoint = "Site uniqueness")

####

alltrends = bind_rows(CPUEtrends, Richtrends, fdistrends, LCBDtrends)

alltrends$Endpoint = as.factor(alltrends$Endpoint)
alltrends$Endpoint = factor(alltrends$Endpoint, levels = levels(alltrends$Endpoint)[c(1,3,4,2)])


alltrends = alltrends %>%
  mutate(Community = ifelse(Type == "Native, non-game",
                            "Local, non-game",
                            ifelse(Type == "Non-native, game",
                                   "Introduced, game",
                            "Whole community")),
         Model = ifelse(grepl("Spatial", Type),
                        "Spatial",
                        "Non-spatial"))

alltrends$Community = as.factor(alltrends$Community)
alltrends$Community = factor(alltrends$Community, levels = levels(alltrends$Community)[c(3,1,2)])

alltrends$Type = as.factor(alltrends$Type)
alltrends$Type = factor(alltrends$Type, levels = levels(alltrends$Type)[c(2,1,3)])

##for the purposes of plotting - specifically adding percent change to the top
##two panels x axis labels; we need to do a lot.
##change labels for the top two panels, bottom two remain the same
##so, I think add 1 and multiple by 10 to get them away from whatever range
##the bottom two are in (bypasses the issue of 0s)

##The odd labeling and breaks and such are all due to the above comment. .plot are 
##the vlaues for plotting purposes only.
alltrends = alltrends %>%
  mutate(Year.trend.plot = ifelse(Endpoint %in% c("Abundance","Rarefied species richness"),
                                                  (Year.trend + 1)*10,
                                                  Year.trend),
         SE.plot = ifelse(Endpoint %in% c("Abundance","Rarefied species richness"),
                          (SE)*10,
                          SE))



breaks_fun <- function(x) {
  if(max(abs(x)) < .026) {
    seq(-.04,.04,.002)
  } else if(max(abs(x)) < .05) {
    seq(-.04,.04,.02)
  } else if(max(abs(x)) > 10.5) {
    seq(9.5,10.5, 0.5)
  } else {seq(9.8,10.2,.2)}
}

limits_fun <- function(x) {
  if (max(abs(x)) < .01) {
    c(-.006,.006)
  } else if(max(abs(x)) < .026) {
    c(-0.028,.028)
  }else if(max(abs(x)) < .05) {
    c(-0.033,0.033)
  }   else if(max(abs(x)) > 10.5) {
    c(9.15,10.85)
  } else {
    c(9.72,10.28)
  }
}

labels_fun <- mylabels <- function(breaks){

  labels <- ifelse(breaks < 1,
         breaks,
         paste(round((breaks/10)-1,2), "\n",ifelse(round(exp((breaks/10)-1),3) > 1,
                                          round(round(exp((breaks/10)-1),3) -1,3),
                                          -(1 - round(exp((breaks/10)-1),3)))*100,"%", sep = ""))

  return(labels)
}

ann_text <- data.frame(Year.trend = c(-.06,.061,
                                      -.017,.018,
                                      -.019,.019,
                                      -.003,.003),
                       Temperature = rep(.55,times = 8),
                       lab = c("Decreasing","Increasing",
                               "Decreasing","Increasing",
                               "More clustering","Less clustering",
                               "Homogenizing","Differentiating"),
                       Endpoint = rep(unique(alltrends$Endpoint), each = 2))

ann_text$Year.trend.plot = ifelse(ann_text$Endpoint %in% c("Abundance",
                                                           "Rarefied species richness"),
                                  (ann_text$Year.trend + 1)*10,
                                  ann_text$Year.trend)
segworkaround = data.frame(Endpoint = rep(unique(alltrends$Endpoint), each = 1),
                           x = c(10,10,0,0),
                           y1 = rep(-Inf, times = 4),
                           y2 = rep(Inf, times = 4))


# %>% filter(!grepl("Spat",Model))
ggplot(alltrends, aes(x = Year.trend.plot, y = Temperature, color = Temperature))+
  facet_wrap(~Endpoint, scales = "free_x")+
  geom_segment(data = segworkaround, inherit.aes = F,
               aes(x = x, xend = x, y = y1,yend = y2), linetype = "dashed")+
  scale_x_continuous(breaks = breaks_fun, limits = limits_fun, labels = labels_fun)+
  ##95%CI
  geom_errorbarh(aes(xmax = Year.trend.plot + SE.plot*1.96, xmin = Year.trend.plot - SE.plot*1.96,
                     group = Type),
                 position = position_dodge(.8),
                 height = 0,
                 linewidth = 2,
                 alpha = 0.5)+
  ##75%CI
  geom_errorbarh(aes(xmax = Year.trend.plot + SE.plot*1.15, xmin = Year.trend.plot - SE.plot*1.15,
                     group = Type),
                 position = position_dodge(.8),
                 height = 0,
                 linewidth = 3,
                 alpha = 0.5)+
  geom_point(aes(
    fill = Community,
    shape = Model,
    group = Type),
    position = position_dodge(.8),
    size = 2,
    # shape = 21,
    color = "black")+
  geom_text(data = ann_text, inherit.aes = F,
            aes(x = Year.trend.plot, y = Temperature, label = lab),
            size = 3)+
  scale_color_manual(values = rev(c("red","violet","blue")), guide= "none")+
  scale_fill_manual(values = rev(c("grey66","white","black")), name = "")+
  scale_alpha_manual(values = c(0.4,1))+
  scale_shape_manual(values = c(21,24), name = "Model")+
  guides(fill = guide_legend(override.aes=list(shape=21),
                             nrow = 1))+
  ylab("Past temperature regime")+
  xlab("Temporal trend per year")+
  theme_bw()+
  theme(axis.text.y = element_text(color = rev(c("red","violet","blue")),
                                   size = 11),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "top",
        legend.box="vertical",
        legend.margin=margin(c(-5,20,-5,-10)),
        axis.text.x = element_text(color = "black",
                                   size = 8),
        axis.title = element_text(size = 12))

ggsave("./Figures/ForestPlotSupp.jpg", dpi = 900, width = 5, height = 5.5, units = "in")


limits_fun <- function(x) {
  if (max(abs(x)) < .01) {
    c(-.005,.005)
  } else if(max(abs(x)) < .026) {
    c(-0.028,.028)
  }else if(max(abs(x)) < .05) {
    c(-0.033,0.033)
  }   else if(max(abs(x)) > 10.5) {
    c(9.15,10.85)
  } else {
    c(9.72,10.28)
  }
}
f2=ggplot(alltrends %>% filter(!grepl("Spat",Model)), aes(x = Year.trend.plot, y = Temperature, color = Temperature))+
  facet_wrap(~Endpoint, scales = "free_x")+
  geom_segment(data = segworkaround, inherit.aes = F,
               aes(x = x, xend = x, y = y1,yend = y2), linetype = "dashed")+
  scale_x_continuous(breaks = breaks_fun, limits = limits_fun, labels = labels_fun)+
  ##95%CI
  geom_errorbarh(aes(xmax = Year.trend.plot + SE.plot*1.96, xmin = Year.trend.plot - SE.plot*1.96,
                     group = Type),
                 position = position_dodge(.6),
                 height = 0,
                 linewidth = 2,
                 alpha = 0.5)+
  ##75%CI
  geom_errorbarh(aes(xmax = Year.trend.plot + SE.plot*1.15, xmin = Year.trend.plot - SE.plot*1.15,
                     group = Type),
                 position = position_dodge(.6),
                 height = 0,
                 linewidth = 3,
                 alpha = 0.5)+
  geom_point(aes(
    fill = Community,
    group = Type),
    position = position_dodge(.6),
    size = 2,
    shape = 21,
    color = "black")+
  geom_text(data = ann_text, inherit.aes = F,
            aes(x = Year.trend.plot, y = Temperature, label = lab),
            size = 3)+
  scale_color_manual(values = rev(c("red","violet","blue")), guide= "none")+
  scale_fill_manual(values = rev(c("grey66","white","black")))+
  guides(fill = guide_legend(override.aes=list(shape=21),
                             nrow = 1))+
  ylab("Past temperature regime")+
  xlab("Temporal trend per year")+
  theme_bw()+
  theme(axis.text.y = element_text(color = rev(c("red","violet","blue")),
                                   size = 11),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "top",
        legend.box="vertical",
        legend.title = element_blank(),
        legend.margin=margin(c(-5,20,-5,-10)),
        axis.text.x = element_text(color = "black",
                                   size = 8),
        axis.title = element_text(size = 12))

ggsave("./Figures/ForestPlotMain_game.jpg",f2, dpi = 900, width = 5.24, height = 5.5, units = "in")
ggsave("./Figures/ForestPlotMain_game.eps",f2, dpi = 900, width = 5.24, height = 5.5, units = "in")

##Figure 4 OPE ----
###trends ----

fOPEdat1 = bind_rows(data.frame(ggemmeans(mOPE,
                                          terms = c("Year[all]","wt_pred_new[13.0,20.5,25.6]"),
                                          weights = "proportional")) %>%
                       mutate(LifeHistory = "Opportunistic"),
                     data.frame(ggemmeans(mOPE2, terms = c("Year[all]","wt_pred_new[13.0,20.5,25.6]"), weights = "proportional")) %>%
                       mutate(LifeHistory = "Periodic"),
                     data.frame(ggemmeans(mOPE3, terms = c("Year[all]","wt_pred_new[13.0,20.5,25.6]"), weights = "proportional")) %>%
                       mutate(LifeHistory = "Equilibrium"))


fOPEdat1$facet = as.factor(fOPEdat1$group )
fOPEdat1$group = as.factor(fOPEdat1$LifeHistory)
fOPEdat1$group = factor(fOPEdat1$group, levels = levels(fOPEdat1$group)[c(1,3,2)])

fOPEdat1$facet2 = ifelse(fOPEdat1$facet == "13",
                         "Cold streams",
                         ifelse(fOPEdat1$facet == "20.5",
                                "Intermediate streams",
                                "Warm streams"))


OPE_trends = ggplot(fOPEdat1, aes(x = x, y = predicted, group = group))+
  facet_grid(~facet2)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = .4,color = NA)+
  geom_line(aes(color = group), size = 0.5)+
  scale_y_continuous(limits = c(.1,.55),
                     expand = c(0,0))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  scale_color_manual(values = rev(c("#d95f02","#59539c","#1b9e77")),
                     labels = rev(c("Opportunistic","Periodic","Equilibrium")),
                     name = "Life history strategy")+
  scale_fill_manual(values = rev(c("#d95f02","#59539c","#1b9e77")),
                    labels = rev(c("Opportunistic","Periodic","Equilibrium")),
                    name = "Life history strategy")+
  ylab("Abundance-weighted assemblage\nlife history continuum")+
  xlab("Year")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 9, color = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        axis.text = element_text(color = "black",
                                 size = 7),
        axis.title = element_text(size = 9),
        legend.margin = margin(t = -5, b = -5))

ggsave("./figures/fig3_lifehistory.png", OPE_trends,dpi = 800, 
       height = 2.5, width = 5)

ggsave("./figures/fig3_lifehistory.eps", OPE_trends,dpi = 800, 
       height = 2.5, width = 5)

##Figure 5 - climate change and introduced game impacts----

###plot----
alldiv = bind_rows(data.frame(ggemmeans(CPUEtrends,
                                        terms = c("estimate[all]","intgame.CPUE.trend[-.02,0,.02]"))) %>%
                     mutate(Div = "Abundance"),
                   data.frame(ggemmeans(Richtrends,
                                        terms = c("estimate[all]","intgame.CPUE.trend[-.02,0,.02]"))) %>%
                     mutate(Div = "Rarefied richness"),
                   data.frame(ggemmeans(LCBDtrends,
                                        terms = c("estimate[all]","intgame.CPUE.trend[-.02,0,.02]"))) %>%
                     mutate(Div = "Site uniqueness"),
                   data.frame(ggemmeans(SEStrends,
                                        terms = c("estimate[all]","intgame.CPUE.trend[-.02,0,.02]"))) %>%
                     mutate(Div = "Functional diversity")
)

alldiv$Div = as.factor(alldiv$Div)
alldiv$Div = factor(alldiv$Div, levels = levels(alldiv$Div)[c(1,3,4,2)])

divcc_ig = ggplot(alldiv, aes(x = x*10, y = predicted*10))+
  facet_wrap(~Div, scales = "free_y")+
  geom_ribbon(aes(ymin = conf.low*10, ymax = conf.high*10, fill = group),
              alpha = 0.25)+
  geom_line(aes(color = group), linewidth = 1)+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  # annotate("text",label = "Increasing\nuniqueness",
  #          y = 0.0005, x = .55, size = 1.5,hjust = 0)+
  # annotate("text",label = "Decreasing\nuniqueness",
  #          y = -0.0004, x = .55, size = 1.5,hjust = 0)+
  # annotate("text",label = "Increasing\ntemperature",
  #          y = -.0044, x = .02, size = 1.5,hjust = 0)+
  # annotate("text",label = "Decreasing\ntemperature",
  #          y = -0.0044, x = -.02, size = 1.5,hjust = 1)+
  scale_x_continuous(limits = c(-.35,.95),
                     breaks = seq(-0.3,.9,.3))+
  scale_color_manual(values = rev(c("#b4e33d","grey63","#fd3e81")),
                     labels = c(-0.2,0,0.2),
                     name = "Introduced, game abundance change per decade")+
  scale_fill_manual(values = rev(c("#b4e33d","grey63","#fd3e81")),
                    labels = c(-0.2,0,0.2),
                    name = "Introduced, game abundance change per decade")+
  ylab("Local, non-game diversity change per decade")+
  xlab(expression(Temperature~change~(degree*C)~per~decade))+
  theme_bw()+
  guides(fill=guide_legend(title.position="top",
                           title.hjust = 0.5),
         color=guide_legend(title.position="top",
                            title.hjust = 0.5))+
  theme(legend.title = element_text(size = 9),
        axis.title = element_text(size = 9),
        axis.text = element_text(color = "black", size = 8),
        strip.text = element_text(size = 9),
        legend.position = "bottom",
        legend.margin = margin(t = -10,
                               b = -5),
        strip.background = element_rect(fill = "white",
                                        color = "black"))


changeplot = ggplot(JulystreamTemps2 %>% filter(estimate*10 > -.5), aes(x = estimate*10, y = PastRegime))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  ggbeeswarm::geom_quasirandom(color = "black", fill = NA,
                               shape = 21,alpha = 0.5,
                               orientation = 'y')+
  geom_boxplot(aes(color = PastRegime), outlier.shape = NA, fill = NA)+
  annotate(geom = "text", x = .7, y = .6,
           label = "Warming", size = 2.5)+
  annotate(geom = "text", x = -.22, y = .6,
           label = "Cooling", size = 2.5)+
  annotate(geom="text", x = -.3, y = 1.25,
           label = "C", size = 3.5)+
  annotate(geom="text", x = -.3, y = 2.25,
           label = "B", size = 3.5)+
  annotate(geom="text", x = -.3, y = 3.25,
           label = "A", size = 3.5)+
  scale_x_continuous(limits = c(-.35,.95),
                     breaks = seq(-0.3,.9,.3))+
  scale_color_manual(values = rev(c("red","violet","blue")))+
  ylab("Past temperature regime")+
  xlab(expression(Temperature~change~(degree*C)~per~decade))+
  # xlab(expression("Temperature change " (degree*C*", 1990-2020"))) +
  theme_bw()+
  theme(axis.text.y = element_text(color = rev(c("red","violet","blue")),
                                   size = 11),
        legend.position = "none",
        axis.text.x = element_text(color = "black",
                                   size = 8),
        axis.title = element_text(size = 9),
        # axis.title.x = element_text(hjust = 1),
        axis.title.y = element_text(hjust = 1),
        plot.margin = margin(t = 0))

plot_grid(plot_grid(NULL, changeplot, NULL,
                    nrow = 1,
                    rel_widths = c(0.3, 1, 0.3),
                    labels = c("","A",""),
                    label_x = -.1),
          # NULL,
          divcc_ig,
          nrow = 2,
          rel_heights = c(.45,1),
          labels = c("","B"))

ggsave("./Figures/Figure5_CC.jpg",
       dpi = 600,
       height = 5.5,
       width = 4,
       units = "in")

ggsave("./Figures/Figure5_CC.eps",
       dpi = 600,
       height = 5.5,
       width = 4,
       units = "in")


 #Supplemental figures----
##figure generation for raw richness and evennes trends is in 4_SupplementalAnalysis.R

##FIgure S2----
fishdatCPUE %>%
  group_by(Year) %>%
  summarize(maxTemp = max(wt_pred_new),
            minTemp = min(wt_pred_new)) %>%
  ungroup() %>%
  summarize(maxTemp = median(maxTemp),
            minTemp = median(minTemp))

temps = ggplot(fishdatCPUE %>%
                 mutate(PastRegime = ifelse(wt_pred_new < lowerbreak,
                                            "Cold",
                                            ifelse(wt_pred_new > upperbreak,
                                                   "Warm",
                                                   "Intermediate")),
                        YearReal = as.character(Year + 1992)),
               aes(x = YearReal, y = wt_pred_new))+
  geom_rect(inherit.aes = F, ymin = 10.8, ymax = 28.5, xmin = -Inf, xmax = Inf, color = NA, fill = "grey90")+
  # geom_violin()+
  ggbeeswarm::geom_quasirandom(aes(color = PastRegime), bandwidth = 0.25, alpha = 0.4)+
  geom_violin(fill = NA)+
  ylab("Mean summer\ntemperature (1990-94)")+
  xlab("Sampling year")+
  scale_x_discrete(breaks = seq(1993,2019,5))+
  # scale_y_continuous(limits = c(0,90))+
  scale_color_manual(name = "Past\ntemperature\nregime",
                     values = c("blue","violet","red"))+
  theme_bw()+
  theme(legend.position = "right",
        legend.margin = margin(c(0,0,0,-5)),
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, color = "black"))


strord = ggplot(fishdatCPUE %>%
                  mutate(StreamOrder = as.numeric(StreamOrder),
                         YearReal = as.character(Year + 1992)) %>%
                  group_by(YearReal, StreamOrder) %>%
                  summarize(count = n()) %>%
                  mutate(count = ifelse(count > 99,
                                        100,
                                        count)),
                aes(x = YearReal, y = StreamOrder))+
  geom_point(aes(size = count, fill = count), shape = 21)+
  ylab("Stream order")+
  xlab("Sampling year")+
  scale_x_discrete(breaks = seq(1993,2019,5))+
  scale_y_continuous(breaks = c(2,4,6,8,10), limits = c(0.25,10))+
  scale_fill_viridis_c(guide = "legend", breaks = c(2,5,25,50,75,100),
                       name = "Sites",
                       labels = c("2","5","25","50","75","100+"))+
  scale_size_continuous(breaks = c(2,5,25,50,75,100), name = "Sites",
                        labels = c("2","5","25","50","75","100+"))+
  # scale_y_continuous(limits = c(0,90))+
  theme_bw()+
  theme(legend.position = "right",
        legend.margin = margin(c(0,0,0,-5)),
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, color = "black"))


landu = ggplot(fishdatCPUE %>%
                 mutate(
                   YearReal = as.character(Year + 1992)) %>%
                 group_by(YearReal, Landuse) %>%
                 summarize(count = n()) %>%
                 mutate(count = ifelse(count > 99,
                                       100,
                                       count)),
               aes(x = YearReal, y = Landuse))+
  geom_point(aes(size = count, fill = count), shape = 21)+
  ylab("Land use")+
  xlab("Sampling year")+
  scale_x_discrete(breaks = seq(1993,2019,5))+
  scale_fill_viridis_c(guide = "legend", breaks = c(2,5,25,50,75,100),
                       name = "Sites",
                       labels = c("2","5","25","50","75","100+"))+
  scale_size_continuous(breaks = c(2,5,25,50,75,100), name = "Sites",
                        labels = c("2","5","25","50","75","100+"))+
  theme_bw()+
  theme(legend.position = "right",
        legend.margin = margin(c(0,0,0,-5)),
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, color = "black"))

fishdatCPUE %>%
  # filter(PredictedWettedWidth_m < 500) %>%
  group_by(Year) %>%
  summarize(maxTemp = max(PredictedWettedWidth_m),
            minTemp = min(PredictedWettedWidth_m)) %>%
  ungroup() %>%
  summarize(maxTemp = median(maxTemp),
            minTemp = median(minTemp))

wetwid = ggplot(fishdatCPUE %>%
                  mutate(YearReal = as.character(Year + 1992)),
                aes(x = YearReal, y = PredictedWettedWidth_m))+
  # geom_violin()+
  geom_rect(inherit.aes = F, ymin = log(2.89), ymax = log(62.3), xmin = -Inf, xmax = Inf, color = NA, fill = "grey90")+
  ggbeeswarm::geom_quasirandom(bandwidth = 0.25, alpha = 0.4)+
  scale_y_continuous(trans = "log", breaks = c(1,5,10,50,100,500,1000))+
  geom_violin(fill = NA)+
  ylab("Predicted wetted\nwidth (m)")+
  xlab("Sampling year")+
  scale_x_discrete(breaks = seq(1993,2019,5))+
  # scale_y_continuous(limits = c(0,90))+
  theme_bw()+
  theme(legend.position = "right",
        legend.margin = margin(c(0,0,0,-5)),
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, color = "black"))

fishdatCPUE %>%
  # filter(WholeConductivity < 10000) %>%
  group_by(Year) %>%
  summarize(maxTemp = max(WholeConductivity),
            minTemp = min(WholeConductivity)) %>%
  ungroup() %>%
  summarize(maxTemp = median(maxTemp),
            minTemp = median(minTemp))

cond = ggplot(fishdatCPUE %>%
                mutate(
                  YearReal = as.character(Year + 1992)),
              aes(x = YearReal, y = WholeConductivity))+
  geom_rect(inherit.aes = F, ymin = log(16.5), ymax = log(1690), xmin = -Inf, xmax = Inf, color = NA, fill = "grey90")+
  ggbeeswarm::geom_quasirandom(bandwidth = 0.25, alpha = .25)+
  scale_y_continuous(trans = "log", breaks = c(1,10,100,1000,10000))+
  geom_violin(fill = NA)+
  ylab(expression(Conductivity~'('*mu*s~cm^'-1'*')'))+
  xlab("Sampling year")+
  scale_x_discrete(breaks = seq(1993,2019,5))+
  # scale_y_continuous(limits = c(0,90))+
  theme_bw()+
  theme(legend.position = "right",
        legend.margin = margin(c(0,0,0,-5)),
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, color = "black"))

yrs = fishdat %>%
  group_by(SiteNumber) %>%
  mutate(minYear = min(Year),
         maxYear = max(Year),
         nyear = maxYear - minYear + 1) %>%
  ungroup() %>%
  mutate(SiteNumber2 = fct_reorder(SiteNumber, minYear)) %>%
  # group_by(SiteNumber)
  mutate(ID = 1:n()) %>%
  # pivot_longer(cols = minYear:maxYear,
  #              names_to = "Type",
  #              values_to = "Year") %>%
  mutate(YearReal = as.character(Year + 1992)) %>%
  ggplot(aes(x = YearReal, y = SiteNumber2), alpha = 0.5)+
  geom_line(aes(group = SiteNumber, color = nyear), linewidth = 0.5)+
  geom_point(shape = 21, fill = NA, alpha = 0.5)+
  scale_color_viridis_c(name = "Years\nsampled")+
  xlab("Sampling year")+
  scale_x_discrete(breaks = seq(1993,2019,5))+
  # scale_y_discrete(expand = c(0,0))+
  theme_bw()+
  # coord_flip()
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "right",
        legend.margin = margin(c(0,0,0,-5)),
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, color = "black"))

strip_x = function(){
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(b = 0,t = 0))
}


cowplot::plot_grid(temps+strip_x(), wetwid+strip_x(),
                   cond+strip_x(), strord+strip_x(),
                   landu+theme(plot.margin = margin(t = 0)),
                   # yrs+theme(plot.margin = margin(t = 0)),
                   ncol = 1,
                   align = "v",
                   axis = "lr",
                   rel_heights = c(.825,.825,.825,.825,1),
                   labels = c(LETTERS[1:5]))

ggsave("./Figures/FigS2.jpg", dpi = 600,
       height = 6.6667,
       width = 6,
       units = "in",
       scale = 1.5)

#Figure S3----
##Non-rarefied Species Richness####
#use glmmTMB, because glmer.nb would not converge
fishdat$SppRichness = rowSums(fishdat[,24:412])

m_rawrich <- glmmTMB(SppRichness ~ Year*HUC2+Agency + StreamOrder + SampleTypeCode +
                       log(PredictedWettedWidth_km) + poly(log(WholeConductivity),2,raw = TRUE) +
                       Year*wt_pred_new+
                       Landuse+Npass+
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
# 13.0  -0.014030 0.00792 Inf  -1.772  0.0765
# 20.5  -0.000675 0.00159 Inf  -0.424  0.6714
# 25.6   0.008405 0.00466 Inf   1.804  0.0712

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
                 Landuse+ Npass+
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
# 13.0  -0.009297 0.00364 Inf  -2.556  0.0106
# 20.5  -0.000696 0.00176 Inf  -0.395  0.6926
# 25.6   0.005153 0.00262 Inf   1.964  0.0495


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

# legendrp = get_plot_component(rawrf1 + theme(legend.box.margin = margin(c(0,0,0,0))),
#                               'guide-box-bottom', return_all = TRUE)
# 
# suppfig2 = plot_grid(plot_grid(rawrf1 + theme(legend.position = "none"), p1,
#                                align = "hv",
#                                nrow = 1,
#                                axis = "l",
#                                labels = c("A","B"),
#                                label_x = .075),
#                      legendrp,
#                      ncol = 1,
#                      rel_heights = c(1,.1))

ggsave("./Figures/FigS2.jpg", suppfig2,dpi = 600, width = 5.75, height = 3, unit = "in")



##trendline figure ----
###Whole Community####
#CPUE
fcdat = data.frame(ggemmeans(m1,
                             terms = c("Year[all]","wt_pred_new[13.0,20.5,25.6]"),
                             weights = "proportional"))

fcdat$group = as.factor(fcdat$group )

CPUE_fig = ggplot(fcdat, aes(x = x, y = predicted, group = group))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_y_continuous(limits = c(0.05,3.1))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  scale_color_manual(values = rev(c("red","violet","blue")),
                     labels = c("Cold","Intermediate","Warm"),
                     name = "Past temperature regime")+
  scale_fill_manual(values = rev(c("red","violet","blue")),
                    labels = c("Cold","Intermediate","Warm"),
                    name = "Past temperature regime")+
  ylab("Abundance (ind./100m min)")+
  xlab("Year")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "bottom",
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9.5),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        plot.margin = margin(b = 2,t=2))

#Rarerich
fr1dat = data.frame(ggemmeans(m_rare,
                              terms = c("Year[all]","wt_pred_new[13.0,20.5,25.6]"),
                              weights = "proportional"))

fr1dat$group = as.factor(fr1dat$group )


rarerich_fig = ggplot(fr1dat, aes(x = x, y = predicted, group = group))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_color_manual(values = rev(c("red","violet","blue")))+
  scale_fill_manual(values = rev(c("red","violet","blue")))+
  scale_y_continuous(limits = c(1.74,8.07))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  ylab("Rarefied species richness")+
  xlab("Year")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9.5),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        plot.margin = margin(b = 2,t=2))

#Functional Diversity

sesdat <- data.frame(ggeffects::ggemmeans(m1f, terms = c("Year[all]","wt_pred_new[13.0,20.5,25.6]"),
                                          weights = "proportional"))
sesdat$group = as.factor(sesdat$group )

ses_fig = ggplot(sesdat, aes(x = x, y = predicted, group = group))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_color_manual(values = rev(c("red","violet","blue")))+
  scale_fill_manual(values = rev(c("red","violet","blue")))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  annotate("text", label = "Overdispersion",
           y = 0.6, x = 6, size = 2.5)+
  annotate("text", label = "Trait clustering",
           y = -1.4, x = 6, size = 2.5)+
  ylab("Functional diversity")+
  xlab("Year")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9.5),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        plot.margin = margin(b = 2,t=2))


##Site uniqueness

lcbddat <- data.frame(ggeffects::ggemmeans(mLCBD, terms = c("Year[all]","wt_pred_new[13.0,20.5,25.6]"),
                                           weights = "proportional"))
lcbddat$group = as.factor(fr1dat$group )

lcbd_fig = ggplot(lcbddat, aes(x = x, y = predicted, group = group))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_y_continuous(limits = c(.2,0.45))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  scale_color_manual(values = rev(c("red","violet","blue")))+
  scale_fill_manual(values = rev(c("red","violet","blue")))+
  ylab("Site uniqueness")+
  annotate("text", label = "More unique composition",
           y = .45, x = 10, size = 2.5)+
  annotate("text", label = "Less unique composition",
           y = .21, x = 10, size = 2.5)+
  xlab("Year")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9.5),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        plot.margin = margin(b = 2,t=2))

#raw rich

frrdat <- data.frame(ggeffects::ggemmeans(m_rawrich, 
                                          terms = c("Year[all]","wt_pred_new[13.0,20.5,25.6]"),
                                          weights = 'proportional'))
frrdat$group = as.factor(frrdat$group )

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
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9.5),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        plot.margin = margin(b = 2,t=2))

#evenness
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
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9.5),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        plot.margin = margin(b = 2,t=2))



###native, non-game####
#CPUE
fr2dat_nat <- data.frame(ggeffects::ggemmeans(m1_native, terms = c("Year[all]","wt_pred_new[13.0,20.5,25.6]"),
                                              weights = "proportional"))
fr2dat_nat$group = as.factor(fr2dat_nat$group )
# fr1dat$group = factor(fr1dat$group, levels= levels(fr1dat$group)[c(2,1,3)] )

CPUE_fig_nng = ggplot(fr2dat_nat, aes(x = x, y = predicted, group = group))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_y_continuous(limits = c(0.05,3.1))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  scale_color_manual(values = rev(c("red","violet","blue")),
                     labels = c("Cold","Intermediate","Warm"),
                     name = "Past temperature regime")+
  scale_fill_manual(values = rev(c("red","violet","blue")),
                    labels = c("Cold","Intermediate","Warm"),
                    name = "Past temperature regime")+
  ylab("Abundance (ind./100m min)")+
  xlab("Year")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9.5),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        plot.margin = margin(b = 2))

#rarefied richness
fr1datnng <- data.frame(ggeffects::ggemmeans(m_rare_native, terms = c("Year[all]","wt_pred_new[13.2, 20.3, 25]"),
                                             weights = "proportional"))
fr1datnng$group = as.factor(fr1dat$group )
# fr1dat$group = factor(fr1dat$group, levels= levels(fr1dat$group)[c(2,1,3)] )

rarerich_fig_nng = ggplot(fr1datnng, aes(x = x, y = predicted, group = group))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_color_manual(values = rev(c("red","violet","blue")))+
  scale_fill_manual(values = rev(c("red","violet","blue")))+
  scale_y_continuous(limits = c(1.74,8.07))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  ylab("Rarefied species richness")+
  xlab("Year")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9.5),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        plot.margin = margin(b = 2))

#SES

sesdat_nng <- data.frame(ggeffects::ggemmeans(m1f_native, terms = c("Year[all]","wt_pred_new[13.0,20.5,25.6]"),
                                              weights = "proportional"))
sesdat_nng$group = as.factor(sesdat_nng$group )

ses_fig_nng = ggplot(sesdat_nng, aes(x = x, y = predicted, group = group))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_color_manual(values = rev(c("red","violet","blue")))+
  scale_fill_manual(values = rev(c("red","violet","blue")))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  annotate("text", label = "Overdispersion",
           y = 0.6, x = 6, size = 2.5)+
  annotate("text", label = "Trait clustering",
           y = -1.4, x = 6, size = 2.5)+
  ylab("Functional diversity")+
  xlab("Year")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9.5),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        plot.margin = margin(b = 2))
#LCBD

lcbddatnng <- data.frame(ggeffects::ggemmeans(mLCBD_native, terms = c("Year[all]","wt_pred_new[13.0,20.5,25.6]"),
                                              weights = "proportional"))
lcbddatnng$group = as.factor(lcbddatnng$group )

lcbd_fig_nng = ggplot(lcbddatnng, aes(x = x, y = predicted, group = group))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_y_continuous(limits = c(.2,0.45))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  scale_color_manual(values = rev(c("red","violet","blue")))+
  scale_fill_manual(values = rev(c("red","violet","blue")))+
  ylab("Site uniqueness")+
  annotate("text", label = "More unique composition",
           y = .45, x = 10, size = 2.5)+
  annotate("text", label = "Less unique composition",
           y = .21, x = 10, size = 2.5)+
  xlab("Year")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9.5),
        plot.margin = margin(b = 2))



###non-native, game####
#CPUE
fr2dat_notnat <- data.frame(ggemmeans(m1_notnative, terms = c("Year[all]","wt_pred_new[13.0,20.5,25.6]"),
                                      weights = "proportional"))
fr2dat_notnat$group = as.factor(fr2dat_notnat$group )
# fr1dat$group = factor(fr1dat$group, levels= levels(fr1dat$group)[c(2,1,3)] )

CPUE_fig_notnat = ggplot(fr2dat_notnat, aes(x = x, y = (predicted), group = group))+
  geom_ribbon(aes(ymin = (conf.low), ymax = (conf.high),fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_y_continuous(limits = c(0.05,3.1))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  scale_color_manual(values = rev(c("red","violet","blue")),
                     labels = c("Cold","Intermediate","Warm"),
                     name = "Past temperature regime")+
  scale_fill_manual(values = rev(c("red","violet","blue")),
                    labels = c("Cold","Intermediate","Warm"),
                    name = "Past temperature regime")+
  ylab("Abundance (ind./100m min)")+
  xlab("Year")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9.5),
        plot.margin = margin(b = 2,r=5))

#rarefied richness
fr1dat_notnat <- data.frame(ggeffects::ggemmeans(m_rare_notnative, terms = c("Year[all]",
                                                                             "wt_pred_new[13.2, 20.3, 25]"),
                                                 weights = "proportional"))
fr1dat_notnat$group = as.factor(fr1dat_notnat$group )
# fr1dat$group = factor(fr1dat$group, levels= levels(fr1dat$group)[c(2,1,3)] )

rarerich_fig_notnat = ggplot(fr1dat_notnat, aes(x = x, y = predicted, group = group))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_color_manual(values = rev(c("red","violet","blue")))+
  scale_fill_manual(values = rev(c("red","violet","blue")))+
  scale_y_continuous(limits = c(1.74,8.07))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  ylab("Rarefied species richness")+
  xlab("Year")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9.5),
        plot.margin = margin(b = 2, r = 5))

#SES

sesdat_notnat <- data.frame(ggeffects::ggemmeans(m1f_notnative, terms = c("Year[all]",
                                                                          "wt_pred_new[13.0,20.5,25.6]"),
                                                 weights = "proportional"))
sesdat_notnat$group = as.factor(sesdat_notnat$group )

ses_fig_notnat = ggplot(sesdat_notnat, aes(x = x, y = predicted, group = group))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_color_manual(values = rev(c("red","violet","blue")))+
  scale_fill_manual(values = rev(c("red","violet","blue")))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  annotate("text", label = "Overdispersion",
           y = 0.6, x = 6, size = 2.5)+
  annotate("text", label = "Trait clustering",
           y = -1.4, x = 6, size = 2.5)+
  ylab("Functional diversity")+
  xlab("Year")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9.5),
        plot.margin = margin(b = 2, r = 5))
#LCBD

lcbddat_notnat <- data.frame(ggeffects::ggemmeans(mLCBD_notnative, terms = c("Year[all]",
                                                                             "wt_pred_new[13.0,20.5,25.6]"),
                                                  weights = "proportional"))
lcbddat_notnat$group = as.factor(lcbddat_notnat$group )

lcbd_fig_notnat = ggplot(lcbddat_notnat, aes(x = x, y = predicted, group = group))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_y_continuous(limits = c(0.2,0.45))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  scale_color_manual(values = rev(c("red","violet","blue")))+
  scale_fill_manual(values = rev(c("red","violet","blue")))+
  ylab("Site uniqueness")+
  annotate("text", label = "More unique composition",
           y = .45, x = 10, size = 2.5)+
  annotate("text", label = "Less unique composition",
           y = .21, x = 10, size = 2.5)+
  xlab("Year")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10, color = "black"),
        legend.position = "none",
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 9.5),
        plot.margin = margin(b = 2, r = 5))


###make final fig####

legend = get_plot_component(CPUE_fig + theme(legend.box.margin = margin(-5,0,-5,0)), 'guide-box-bottom', return_all = TRUE)
CPUE_fig2 <- CPUE_fig + theme(legend.position='none')
# CPUE_fig2, rarerich_fig, ses_fig, lcbd_fig,
# hist1,

##remake this figure:
#3 columns: whole, native, non-native
#6 rows: CPUE, Rarefich, LCBD, funcdiv (for all); raw rich and evenness for whole
##native and non-native can all get rid of axis title and text for y (leave ticks)
##funcdiv for native and non-native, evenness for whole should have x axes
##  all others can get rid of x-axis text and titles
##Lettering goes reading format (right then down)
##ABC
##DEF
##GHI
##JKL
##M00
##N00
##legend

theme_noxy = function(){
  theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
}
theme_nox = function(){
  theme(axis.title.x = element_blank()
        )
}

figunivar_nng <- plot_grid(CPUE_fig2 + theme_nox(), CPUE_fig_nng + theme_noxy(), CPUE_fig_notnat+ theme_noxy(),
                           rarerich_fig+ theme_nox(), rarerich_fig_nng+ theme_noxy(), rarerich_fig_notnat+ theme_noxy(),
                           lcbd_fig+ theme_nox(), lcbd_fig_nng+ theme_noxy(), lcbd_fig_notnat+ theme_noxy(),
                           ses_fig+ theme_nox(), ses_fig_nng+ theme_noxy(), ses_fig_notnat+ theme_noxy(),
                           rawrf1+ theme_nox(), NULL, NULL,
                           p1+ theme_nox(), NULL, NULL,
                           labels = c(LETTERS[1:13],"","","N","",""),
                           label_size = 15,
                           label_x = c(-.025,-.025,-.015,
                                       -.025,-.025,-.015,
                                       -.025,-.025,-.015,
                                       -.025,-.025,-.015,
                                       .05,0,0,
                                       -.025,0,0),
                           label_y = c(1,1,1,
                                       1,1,1,
                                       1,1,1,
                                       1,1,1,
                                       1.05,1,1,
                                       1,1,1),
                           align = "hv",
                           axis = "lr",
                           nrow = 6,
                           ncol = 3)



figunivar_nng_l = plot_grid(figunivar_nng, legend, ncol = 1, rel_heights = c(1, .1))

ggsave("./figures/figUnivar_comparisons_YEARC_LUincluded.png", figunivar_nng_l, dpi = 900, 
       height = 8.5, width = 4.5, scale = 1.5)



# figunivar_nng <- plot_grid(CPUE_fig2, rarerich_fig, lcbd_fig, ses_fig,
#                            CPUE_fig_nng, rarerich_fig_nng,lcbd_fig_nng, ses_fig_nng,
#                            CPUE_fig_notnat, rarerich_fig_notnat,lcbd_fig_notnat, ses_fig_notnat,
#                            labels = c("A", "D", "G", "J",
#                                       "B", "E", "H", "K",
#                                       "C", "F", "I", "L"),
#                            label_size = 15,
#                            label_x = c(-.025,-.025,-.025,0,
#                                        -.025,-.025,-.025,0,
#                                        -.025,-.025,-.015,0),
#                            rel_heights = c(.825,.825,1),
#                            align = "v",
#                            axis = "lr",
#                            nrow = 3,
#                            ncol = 4)





####Landuse






##landuse figure----

#CPUE
fcdatLU = data.frame(ggemmeans(m1, terms = c("Year[all]","Landuse"),
                               weights = "proportional"))

fcdatLU$group = as.factor(fcdatLU$group )
fcdatLU$group = factor(fcdatLU$group, levels = levels(fcdatLU$group)[c(1,3,4,2)])
levels(fcdatLU$group)[1] = "Forest/wetland"
levels(fcdatLU$group)[2] = "Grassland/shrub"


CPUE_figLU = ggplot(fcdatLU, aes(x = x, y = predicted, group = group))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  geom_line(aes(color = group),size = 1)+
  scale_y_continuous(limits = c(0.05,3.02))+
  scale_x_continuous(breaks = c(1,10,20,27),
                     labels = c(1993,2003,2013,2019))+
  scale_color_manual(values = (c("#1A9B75","#C93D80","#d85e02","#746DB7")))+
  scale_fill_manual(values = (c("#1A9B75","#C93D80","#d85e02","#746DB7")))+
  ylab("Abundance (ind./100m min)")+
  xlab("Year")+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black",
                             size = 8),
    axis.title = element_text(size = 9))


#richness
fr1datLU = data.frame(ggemmeans(m_rare, terms = c("Landuse"),
                                weights = "proportional"))


fr1datLU$x = as.factor(fr1datLU$x)
fr1datLU$x = factor(fr1datLU$x, levels = levels(fr1datLU$x)[c(1,3,4,2)])
levels(fr1datLU$x)[1] = "Forest/wetland"
levels(fr1datLU$x)[2] = "Grassland/shrub"

labs_lurich = data.frame(x = fr1datLU$x,
                         y = fr1datLU$conf.high+.1,
                         labs = c("A","B","A","A"))

rarerich_figLU = ggplot(fr1datLU, aes(x = x, y = predicted, color = x))+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.25)+
  geom_point()+
  geom_text(data = labs_lurich, inherit.aes = F,
            aes(x = x, y = y, label = labs),
            size = 3)+
  # geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  # geom_line(aes(color = group),size = 1)+
  # scale_y_continuous(limits = c(0.1,4.5))+
  # scale_x_continuous(breaks = c(1,10,20,27),
  #                    labels = c(1993,2003,2013,2019))+
  scale_color_manual(values = (c("#1A9B75","#C93D80","#d85e02","#746DB7")))+
  scale_fill_manual(values = (c("#1A9B75","#C93D80","#d85e02","#746DB7")))+
  ylab("Rarefied species richness")+
  xlab("")+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black",
                             size = 8),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank())

#SES

sesdatLU = data.frame(ggemmeans(m1f, terms = c("Landuse"),
                                weights = "proportional"))


sesdatLU$x = as.factor(sesdatLU$x)
sesdatLU$x = factor(sesdatLU$x, levels = levels(sesdatLU$x)[c(1,3,4,2)])
levels(sesdatLU$x)[1] = "Forest/wetland"
levels(sesdatLU$x)[2] = "Grassland/shrub"

#forest and grass different from ag, not other diffs, so
#forest: B, Grass: B, Ag: A, Urb: AB

###CHECK THAT THIS MATCHES
labs_luses = data.frame(x = sesdatLU$x,
                         y = sesdatLU$conf.high+.025,
                         labs = c("B","AB","B","A"))

ses_figLU = ggplot(sesdatLU, aes(x = x, y = predicted, color = x))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.25)+
  geom_point()+
  geom_text(data = labs_luses, inherit.aes = F,
            aes(x = x, y = y, label = labs),
            size = 3)+
  # geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  # geom_line(aes(color = group),size = 1)+
  # scale_y_continuous(limits = c(0.1,4.5))+
  # scale_x_continuous(breaks = c(1,10,20,27),
  #                    labels = c(1993,2003,2013,2019))+
  scale_color_manual(values = (c("#1A9B75","#C93D80","#d85e02","#746DB7")))+
  scale_fill_manual(values = (c("#1A9B75","#C93D80","#d85e02","#746DB7")))+
  ylab("Functional diversity")+
  xlab("")+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black",
                             size = 8),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank())


#lcbd

lcbddatLU = data.frame(ggemmeans(mLCBD, terms = c("Landuse"),
                                 weights = "proportional"))

lcbddatLU$x = as.factor(lcbddatLU$x)
lcbddatLU$x = factor(lcbddatLU$x, levels = levels(lcbddatLU$x)[c(1,3,4,2)])
levels(lcbddatLU$x)[1] = "Forest/wetland"
levels(lcbddatLU$x)[2] = "Grassland/shrub"

#Forest and grass different from Ag, no other diffs,
#forest: B, Grass: B, Ag: A, Urb: AB
labs_lulcbd = data.frame(x = lcbddatLU$x,
                         y = lcbddatLU$conf.high+.005,
                         labs = c("B","AB","B","A"))

lcbd_figLU = ggplot(lcbddatLU, aes(x = x, y = predicted, color = x))+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.25)+
  geom_point()+
  geom_text(data = labs_lulcbd, inherit.aes = F,
            aes(x = x, y = y, label = labs),
            size = 3)+
  # geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill = group),alpha = 0.25,color = NA)+
  # geom_line(aes(color = group),size = 1)+
  # scale_y_continuous(limits = c(0.1,4.5))+
  # scale_x_continuous(breaks = c(1,10,20,27),
  #                    labels = c(1993,2003,2013,2019))+
  scale_color_manual(values = (c("#1A9B75","#C93D80","#d85e02","#746DB7")))+
  scale_fill_manual(values = (c("#1A9B75","#C93D80","#d85e02","#746DB7")))+
  ylab("Site uniqueness")+
  xlab("")+
  theme_bw()+
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black",
                             size = 8),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title.x = element_blank())


figunivarLU <- cowplot::plot_grid(CPUE_figLU, rarerich_figLU, lcbd_figLU, ses_figLU,
                                  labels = c("A","B", "C","D"),
                                  label_size = 12.5,
                                  align = "h",
                                  axis = "b",
                                  nrow = 1,
                                  ncol = 4,
                                  label_y = 1)

ggsave("./figures/figUnivar_Landuse.jpg", figunivarLU, dpi = 800, 
       height = 1.65, width = 6, scale = 1.5)

##OPE fish characterists----
###local vs introduced,game ----

OPE_n = bind_rows(data.frame(ggemmeans(mO, terms =  c("Native"), weights = "proportional")) %>%
                    mutate(LifeHistory = "Opportunistic"),
                  data.frame(ggemmeans(mE, terms =  c("Native"), weights = "proportional")) %>%
                    mutate(LifeHistory = "Equilibrium"),
                  data.frame(ggemmeans(mP, terms =  c("Native"), weights = "proportional")) %>%
                    mutate(LifeHistory = "Periodic"))

OPENvNN = ggplot(OPE_n, aes(x  = x, y = predicted, color = LifeHistory))+
  # facet_wrap(~group)+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.25))+
  geom_point(size = 2.5,
             position = position_dodge(width = 0.25))+
  ylab("Life history continuum")+
  scale_y_continuous(limits = c(.1,.65),
                     expand = c(0,0))+
  scale_x_discrete(labels = c("Local",
                              "Introduced"))+
  scale_color_manual(values = (c("#1b9e77","#d95f02","#59539c")),
                     labels = (c("Equilibrium","Opportunistic", "Periodic")),
                     name = "Life history strategy")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black",
                                 size = 7),
        axis.title = element_text(size = 9))


OPE_g = bind_rows(data.frame(ggemmeans(mO, terms =  c("Status"), weights = "proportional")) %>%
                    mutate(LifeHistory = "Opportunistic"),
                  data.frame(ggemmeans(mE, terms =  c("Status"), weights = "proportional")) %>%
                    mutate(LifeHistory = "Equilibrium"),
                  data.frame(ggemmeans(mP, terms =  c("Status"), weights = "proportional")) %>%
                    mutate(LifeHistory = "Periodic"))

OPEGvNG = ggplot(OPE_g, aes(x  = x, y = predicted, color = LifeHistory))+
  # facet_wrap(~group)+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0,
                position = position_dodge(width = 0.25))+
  geom_point(size = 2.5,
             position = position_dodge(width = 0.25))+
  ylab("Life history continuum")+
  scale_y_continuous(limits = c(.1,.65),
                     expand = c(0,0))+
  # scale_x_discrete(labels = c("Native,\nnon-game",
  #                             "Non-native,\ngame"))+
  scale_color_manual(values = (c("#1b9e77","#d95f02","#59539c")),
                     labels = (c("Equilibrium","Opportunistic", "Periodic")),
                     name = "Life history strategy")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(color = "black",
                                 size = 7),
        axis.title = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = margin(l = 5,r=1,b=5,t=5))


###Opportunistic/Periodic and temp####

fishstats_t = fishstats %>%
  left_join(FishLifeTraits %>%
              dplyr::select(Taxa, temperature) %>%
              filter(Taxa %in% fishstats$species),
            by = join_by(species == Taxa))


fishstats_tl = fishstats_t %>%
  pivot_longer(cols = Equilibrium:Opportunist,
               names_to = "LifeHistory",
               values_to = "Continuum") %>%
  filter(LifeHistory != "Equilibrium")

OPE_temp = ggplot(fishstats_tl, aes(y = temperature, x = Continuum))+
  facet_wrap(~LifeHistory)+
  # geom_boxplot()+
  geom_point(shape = 21, fill = NA, size = 2)+
  stat_smooth(aes(color = LifeHistory, fill = LifeHistory), method = "lm")+
  scale_color_manual(values = c("#d95f02","#59539c"),
                     labels = (c("Opportunistic", "Periodic")))+
  scale_fill_manual(values = c("#d95f02","#59539c"),
                    labels = (c("Opportunistic", "Periodic")))+
  xlab("Life history continuum")+
  ylab("Fish temperature preference")+
  scale_y_continuous(labels = scales::label_number(suffix = "°C"),
                     limits = c(5,26))+
  theme_bw()+
  theme(        axis.text = element_text(color = "black",
                                         size = 7),
                axis.title = element_text(size = 9),
                axis.title.y = element_text(hjust = -.2),
                legend.position = "none",
                strip.background = element_blank(),
                strip.text = element_blank()) 


###make full figure ####

legendA = get_plot_component(OPEGvNG + theme(legend.box.margin = margin(c(0,0,0,0))), 'guide-box-bottom', return_all = TRUE)


bottom_row = plot_grid(OPENvNN, 
                       OPEGvNG+ theme(legend.position = "none"),
                       OPE_temp,
                       align = "h",
                       axis = "b",
                       nrow = 1,
                       labels = c("A","B","C"),
                       rel_widths = c(.55,.425,1),
                       label_x = c(0,-.1,.05),
                       label_y = 1.025)

fullsuppfig6 = plot_grid(bottom_row, legendA,
          nrow = 2,
          rel_heights = c(1,.1))

ggsave("./figures/supfig_lifehistory.png", fullsuppfig6,dpi = 800, 
       height = 2.5, width = 6.5)



##site-level trend map####
##for each endpoint, generate trends for sequential temperature, min:max, by = 0.01
##then match to dataset by rounding wt_pred_new to nearest 0.01
##then map by endpoint trend, facet_wrap(~endpoint)
##
wt_predrange = seq(round(min(fishdatCPUE$wt_pred_new),1), round(max(fishdatCPUE$wt_pred_new),1),by = 0.1)

##CPUE needs to be based on landuse as well, since landuse interacts with time for this model
CPUEtrends = emtrends(m1, ~wt_pred_new|Landuse, var = "Year",
         at = list(wt_pred_new = wt_predrange),lmer.df = "asymptotic",
         weights = "proportional"
         )
richnesstrends = emtrends(m_rare, ~wt_pred_new, var = "Year",
                      at = list(wt_pred_new = wt_predrange),lmer.df = "asymptotic",
                      weights = "proportional"
)
SEStrends = emtrends(m1f, ~wt_pred_new, var = "Year",
                      at = list(wt_pred_new = wt_predrange),lmer.df = "asymptotic",
                     weights = "proportional"
)
LCBDtrends = emtrends(mLCBD, ~wt_pred_new, var = "Year",
                      at = list(wt_pred_new = wt_predrange),lmer.df = "asymptotic",
                      weights = "proportional"
)

trends = data.frame(CPUEtrends) %>%
  rename(CPUE.trend = Year.trend) %>%
  dplyr::select(Landuse, wt_pred_new,  CPUE.trend) %>%
  left_join(data.frame(richnesstrends) %>%
              rename(richness.trend = Year.trend) %>%
              dplyr::select(wt_pred_new,  richness.trend)) %>%
  left_join(data.frame(SEStrends) %>%
              rename(SES.trend = Year.trend) %>%
              dplyr::select(wt_pred_new,SES.trend)) %>%
  left_join(data.frame(LCBDtrends) %>%
              rename(LCBD.trend = Year.trend) %>%
              dplyr::select(wt_pred_new,  LCBD.trend)) %>%
  mutate(wt_pred_new = as.character(wt_pred_new))

temporaltrendwsite = fishdat %>%
  group_by(SiteNumber) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(wt_pred_new = round(wt_pred_new, 1))  %>%
  mutate(wt_pred_new = as.character(wt_pred_new)) %>%
  dplyr::select(Latitude_dd, Longitude_dd, wt_pred_new, Landuse) %>%
  left_join(trends)
  
temporaltrendwsite = temporaltrendwsite %>%
  pivot_longer(cols = contains("trend"),
               names_to = "endpoint",
               values_to = "trend") %>%
  mutate(endpoint = gsub(".trend","",endpoint),
         endpoint = ifelse(endpoint == "richness",
                           "Rarefied species richness",
                           ifelse(endpoint == "CPUE",
                                  "Abundance",
                                  ifelse(endpoint == "LCBD",
                                         "Site uniqueness",
                                         "Functional diversity"))))

##started with the pts dataset, instead of temporaltrendwsite
world1 <- sf::st_as_sf(maps::map('state', plot = FALSE, fill = TRUE))

pts <- readRDS("./Data/sitesLatLong 3.rds")
pts <- pts[!is.na(pts$wt_pred_new),]
pts <- st_as_sf(pts, coords = c("Longitude_dd","Latitude_dd"), crs = 4269)
pts <- st_transform(pts,crs = 5070)

pts = pts %>%
  mutate(wt_pred_new = as.character(round(wt_pred_new,1))) %>%
  left_join(trends %>% mutate(wt_pred_new = as.character(wt_pred_new))) %>%
  mutate(wt_pred_new = as.numeric(wt_pred_new)) %>%
  pivot_longer(cols = contains("trend"),
               names_to = "endpoint",
               values_to = "Estimate") %>%
  mutate(endpoint = gsub(".trend","",endpoint),
         endpoint = ifelse(endpoint == "richness",
                           "Rarefied species richness",
                           ifelse(endpoint == "CPUE",
                                  "Abundance",
                                  ifelse(endpoint == "LCBD",
                                         "Site uniqueness",
                                         "Functional diversity"))))
pts %>% arrange(endpoint)

endpoint_mapshld = pts %>%
  split(.$endpoint) %>%
  map(~ ggplot(.)+
        facet_wrap(~endpoint)+
        geom_sf(data = world1)+
        geom_sf(aes(fill = Estimate), shape = 21, size=1, alpha = 0.6)+
        scale_fill_gradient2(midpoint = 0,
                             high = ("#b4e33d"),
                             low = ("#fd3e81"),
                             mid = "white")+
        theme_bw()+
        theme(legend.position = "none",
              strip.text = element_text(size = 9),
              axis.text = element_text(size = 7),
              plot.margin = margin(0,-1,0,-1)))

endpoint_maps = cowplot::plot_grid(endpoint_mapshld$Abundance,
                                   endpoint_mapshld$`Rarefied species richness`,
                                   endpoint_mapshld$`Site uniqueness`,
                                   endpoint_mapshld$`Functional diversity`)

lgend = pts %>%
  filter(endpoint == "Site uniqueness") %>%
  ggplot()+
        # facet_wrap(~endpoint)+
        geom_sf(data = world1)+
        geom_sf(aes(fill = Estimate), shape = 21, size=1, alpha = 0.6)+
        scale_fill_gradient2(midpoint = 0,
                             high = ("#b4e33d"),
                             low = ("#fd3e81"),
                             mid = "white",
                             limits = c(-0.004,0.004),
                             name = "",
                             breaks = c(-0.004,0,0.004),
                             labels = c("Decreasing diversity",
                                        "No change", "Increasing diversity"))+
        theme_bw()+
        theme(legend.position = "bottom",
              legend.key.height = unit(0.15, 'cm'),
              legend.key.width = unit(1.5, "cm"),
              legend.title = element_blank()
              )

legendA = get_plot_component(lgend + theme(legend.box.margin = margin(c(0,0,0,0))), 'guide-box-bottom', return_all = TRUE)
      
      
endpoint_maps2 = plot_grid(endpoint_maps,
                          legendA, nrow = 2, rel_heights = c(1,0.15))

ggsave("./Figures/EndpointMap.jpg", endpoint_maps2,dpi = 900,
       height = 4.5,
       width = 6,
       units = "in")



##temporal figure of past temperature regime histogram ####

ggplot(fishdat %>%
         mutate(PastRegime = ifelse(wt_pred_new < lowerbreak,
                                    "Cold",
                                    ifelse(wt_pred_new > upperbreak,
                                           "Warm",
                                           "Intermediate")),
                YearReal = Year + 1992),
       aes(x = wt_pred_new, fill = PastRegime))+
  facet_wrap(~YearReal)+
  geom_histogram(binwidth = 1)+
  xlab("Mean summer temperature (1990-94)")+
  ylab("Number of sites")+
  scale_y_continuous(limits = c(0,90))+
  scale_fill_manual(name = "Past temperature regime",
                    values = c("blue","violet","red"))+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.margin = margin(c(-5,0,0,0)),
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, color = "black"))

ggsave("./Figures/TempHistogramThruTime.jpg", dpi = 300,
       width = 6,
       height = 5.5,
       units = "in")


##temporal figure of stream order histogram ####


ggplot(fishdat %>% mutate(StreamOrder = as.numeric(StreamOrder)),aes(x = StreamOrder))+
  facet_wrap(~CollectionYear)+
  geom_histogram(binwidth = 1, fill = "grey78", color = "black")+
  xlab("Stream order")+
  ylab("Number of sites")+
  scale_y_continuous(limits = c(0,165))+
  scale_x_continuous(breaks = 1:10)+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.margin = margin(c(-5,0,0,0)),
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8, color = "black"))

ggsave("./Figures/StreamOrderHistorgrtmapamThruTime.jpg", dpi = 300,
       width = 6,
       height = 5.5,
       units = "in")