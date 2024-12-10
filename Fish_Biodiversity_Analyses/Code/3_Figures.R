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
               text.size = .7, lwd = 1, width = .25)+
  tm_layout(frame = FALSE, 
            inner.margins = c(0.04,0.0,0.0,0),
            asp = 1.6,
            bg.color = "transparent")

##make sure plot window is sized so that ratio is 1.1
MAPHOLD = tmap_grob(Plot_Thermal)

fullfig1 = plot_grid(plot_grid(NULL,hist12,NULL,
                               nrow = 3,
                               rel_heights = c(0.25,1,0.25),
                               labels = c("","A",""),
                               label_y = 1.1),
                     # NULL,
                     MAPHOLD,
          nrow = 1,
          labels = c("","B"),
          rel_widths = c(.8,1),
          label_x = 0,
          label_y = .9)

ggsave("./Figures/Figure1AB_new.jpg",fullfig1,
       height = 3, width = 5, units = "in",
       dpi = 600)


##Figure 3 Past stream temp forest plot######
#CPUE
test(emtrends(m1, ~wt_pred_new, var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE    df t.ratio p.value
# 13.0   -0.02325 0.00941 455.9  -2.471  0.0138
# 20.5    0.00483 0.00502  42.6   0.962  0.3413
# 25.6    0.02393 0.00717 183.1   3.336  0.0010

test(emtrends(m1_native,  ~ wt_pred_new,
              var = "Year", at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE    df t.ratio p.value
# 13.0   -0.05500 0.01249 597.4  -4.404  <.0001
# 20.5    0.00202 0.00624  47.3   0.324  0.7476
# 25.6    0.04079 0.00939 255.4   4.343  <.0001

test(emtrends(m1.sp,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE  df z.ratio p.value
# 13.0   -0.02022 0.00918 Inf  -2.203  0.0276
# 20.5    0.00249 0.00473 Inf   0.526  0.5988
# 25.6    0.01793 0.00693 Inf   2.588  0.0097

test(emtrends(m1_notnative,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE  df t.ratio p.value
# 13.0    0.00952 0.01039 347   0.916  0.3603
# 20.5    0.01215 0.00584  39   2.082  0.0440
# 25.6    0.01394 0.00775 125   1.799  0.0744

CPUEtrends = data.frame(Year.trend = c(-0.02325,0.00483,0.02393,
                                       -0.05500,0.00202,0.04079,
                                       -0.02022,-0.00249,0.01793,
                                       0.00952,0.01215,0.01394),
                        SE = c(0.00941,0.00502,0.00717,
                               0.01249,0.00624,0.00939,
                               0.00918,0.00473,0.00693,
                               0.01039,0.00584,0.00775),
                        Temperature = rep(c("Cold","Intermediate","Warm"), times = 4),
                        Type = rep(c("Whole community","Native, non-game","Spatial model","Non-native, game"), each = 3),
                        Endpoint = "Abundance")

CPUEtrends$decadepercchange = exp(CPUEtrends$Year.trend)^10 - 1


#richness
test(emtrends(m_rare, ~wt_pred_new, var = "Year", pbkrtest.limit = 60000,
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE    df t.ratio p.value
# 13.0   -0.01430 0.00288 880.5  -4.974  <.0001
# 20.5   -0.00261 0.00121  42.5  -2.151  0.0372
# 25.6    0.00534 0.00196 322.1   2.726  0.0068

test(emtrends(m_rare_native, ~wt_pred_new, var = "Year", pbkrtest.limit = 60000,
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE     df t.ratio p.value
# 13.0   -0.00924 0.00383 1485.8  -2.414  0.0159
# 20.5   -0.00305 0.00124   33.7  -2.452  0.0195
# 25.6    0.00115 0.00220  358.4   0.524  0.6005


test(emtrends(m_rare.sp,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE  df z.ratio p.value
# 13.0   -0.01369 0.00282 Inf  -4.862  <.0001
# 20.5   -0.00277 0.00109 Inf  -2.542  0.0110
# 25.6    0.00465 0.00189 Inf   2.460  0.0139

test(emtrends(m_rare_notnative, ~wt_pred_new, var = "Year", pbkrtest.limit = 60000,
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE    df t.ratio p.value
# 13.0   -0.01391 0.00376 653.3  -3.695  0.0002
# 20.5   -0.00424 0.00181  45.8  -2.346  0.0234
# 25.6    0.00234 0.00256 185.5   0.912  0.3628

Richtrends = data.frame(Year.trend = c(-0.01430,-0.00261,0.00534,
                                       -0.00924,-0.00305,0.00115,
                                       -0.01369,-0.00277,0.00465,
                                       -0.01391,-0.00424,0.00234),
                        SE = c(0.00288,0.00121,0.00196,
                               0.00383,0.00124,0.00220,
                               0.00282,0.00109,0.00189,
                               0.00376,0.00181,0.00256),
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
# 13.0    0.00886 0.00800 871.2   1.107  0.2685
# 20.5    0.00246 0.00349  40.8   0.705  0.4850
# 25.6   -0.00189 0.00529 234.6  -0.357  0.7214

test(emtrends(m1f_native,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE     df t.ratio p.value
# 13.0  -0.000292 0.00748 1006.8  -0.039  0.9689
# 20.5  -0.002539 0.00296   34.8  -0.858  0.3966
# 25.6  -0.004067 0.00484  300.2  -0.841  0.4010

test(emtrends(m1f.sp,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE  df z.ratio p.value
# 13.0    0.01460 0.00738 Inf   1.979  0.0478
# 20.5    0.00191 0.00268 Inf   0.715  0.4744
# 25.6   -0.00671 0.00490 Inf  -1.370  0.1706

test(emtrends(m1f_notnative, ~wt_pred_new,
              var = "Year", 
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend      SE    df t.ratio p.value
# 13.0    0.01413 0.00815 935.2   1.735  0.0830
# 20.5   -0.00119 0.00337  35.9  -0.353  0.7261
# 25.6   -0.01161 0.00530 256.9  -2.189  0.0295

fdistrends = data.frame(Year.trend = c(0.00886,0.00246,-0.00189,
                                       -0.000292,-0.002539,-0.004067,
                                       0.01460,0.00191,-0.00671,
                                       0.01413,-0.00119, -0.01161),
                        SE = c(0.00800,0.00349,0.00529,
                               0.00748,0.00296,0.00484,
                               0.00738,0.00268,0.00490,
                               0.00815,0.00337,0.00530),
                        Temperature = rep(c("Cold","Intermediate","Warm"), times = 4),
                        Type = rep(c("Whole community","Native, non-game","Spatial model","Non-native, game"), each = 3),
                        Endpoint = "Functional diversity")

#LCBD

test(emtrends(mLCBD,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend       SE    df t.ratio p.value
# 13.0   0.003020 0.000816 531.8   3.702  0.0002
# 20.5  -0.000615 0.000408  43.9  -1.507  0.1389
# 25.6  -0.003087 0.000587 201.6  -5.262  <.0001

test(emtrends(mLCBD_native,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend       SE     df t.ratio p.value
# 13.0   0.000480 0.001150 1005.6   0.417  0.6766
# 20.5  -0.000524 0.000488   48.7  -1.075  0.2877
# 25.6  -0.001207 0.000715  242.3  -1.688  0.0927

test(emtrends(mLCBD.sp,  ~ wt_pred_new,
              var = "Year",
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend       SE  df z.ratio p.value
# 13.0   0.002634 0.000759 Inf   3.471  0.0005
# 20.5  -0.000292 0.000346 Inf  -0.843  0.3994
# 25.6  -0.002282 0.000539 Inf  -4.231  <.0001

test(emtrends(mLCBD_notnative, ~wt_pred_new, var = "Year", pbkrtest.limit = 60000,
              at = list(wt_pred_new = c(13.0,20.5,25.6)),
              weights = "proportional"))
# wt_pred_new Year.trend       SE    df t.ratio p.value
# 13.0    0.00259 0.000944 672.5   2.742  0.0063
# 20.5   -0.00116 0.000400  37.7  -2.904  0.0061
# 25.6   -0.00371 0.000646 354.7  -5.741  <.0001

LCBDtrends = data.frame(Year.trend = c(0.003020,-0.000615,-0.003087,
                                       0.000480,-0.000524,-0.001207,
                                       0.002634,-0.000292,-0.002282,
                                       0.00259,-0.00116,-0.00371),
                        SE = c(0.000816,0.000408,0.000587,
                               0.001150,0.000488,0.000715,
                               0.000759,0.000346,0.000539,
                               0.000944,0.000400,0.000646
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
  xlab("Temporal trend")+
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
  xlab("Temporal trend")+
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
  ylab("Assemblage life history continuum")+
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
  ylab("Diversity change per decade")+
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

 #Supplemental figures----
##figure generation for raw richness and evennes trends is in 4_SupplementalAnalysis.R

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
  scale_y_continuous(limits = c(0.05,3))+
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
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
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
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
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
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
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
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
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
  scale_y_continuous(limits = c(0.05,3))+
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
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
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
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
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
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
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
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
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
  scale_y_continuous(limits = c(0.05,3))+
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
        plot.margin = margin(b = 2))

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
        plot.margin = margin(b = 2))

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
        plot.margin = margin(b = 2))
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
        plot.margin = margin(b = 2))


###make final fig####

legend = get_plot_component(CPUE_fig + theme(legend.box.margin = margin(-5,0,-5,0)), 'guide-box-bottom', return_all = TRUE)
CPUE_fig2 <- CPUE_fig + theme(legend.position='none')
# CPUE_fig2, rarerich_fig, ses_fig, lcbd_fig,
# hist1,

figunivar_nng <- plot_grid(CPUE_fig2, rarerich_fig, lcbd_fig, ses_fig,
                           CPUE_fig_nng, rarerich_fig_nng,lcbd_fig_nng, ses_fig_nng,
                           CPUE_fig_notnat, rarerich_fig_notnat,lcbd_fig_notnat, ses_fig_notnat,
                           labels = c("A", "D", "G", "J",
                                      "B", "E", "H", "K",
                                      "C", "F", "I", "L"),
                           label_size = 15,
                           label_x = c(-.025,-.025,-.025,0,
                                       -.025,-.025,-.025,0,
                                       -.025,-.025,-.015,0),
                           rel_heights = c(.825,.825,1),
                           align = "v",
                           axis = "lr",
                           nrow = 3,
                           ncol = 4)

figunivar_nng_l = plot_grid(figunivar_nng, legend, ncol = 1, rel_heights = c(1, .1))

ggsave("./figures/figUnivar_comparisons_YEARC_LUincluded.png", figunivar_nng_l, dpi = 900, 
       height = 4.25, width = 6, scale = 1.5)




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
  scale_y_continuous(limits = c(0.05,3))+
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
        legend.position = "none",
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

bottom_row = plot_grid(OPENvNN, OPEGvNG,OPE_temp,
                       align = "h",
                       axis = "b",
                       nrow = 1,
                       labels = c("A","B","C"),
                       rel_widths = c(.55,.425,1),
                       label_x = c(0,-.1,.05),
                       label_y = 1.025)

ggsave("./figures/supfig_lifehistory.png", bottom_row,dpi = 800, 
       height = 2, width = 6.5)



##site-level trend map####
##for each endpoint, generate trends for sequential temperature, min:max, by = 0.01
##then match to dataset by rounding wt_pred_new to nearest 0.01
##then map by endpoint trend, facet_wrap(~endpoint)
##
wt_predrange = seq(round(min(fishdatCPUE$wt_pred_new),1), round(max(fishdatCPUE$wt_pred_new),1),by = 0.1)

CPUEtrends = emtrends(m1, ~wt_pred_new, var = "Year",
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
  dplyr::select(wt_pred_new,  CPUE.trend) %>%
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
  dplyr::select(Latitude_dd, Longitude_dd, wt_pred_new) %>%
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
