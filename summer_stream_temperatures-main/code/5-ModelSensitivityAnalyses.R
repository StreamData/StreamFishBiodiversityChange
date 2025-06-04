final_splm = read_rds('./data/splm_selected.2024.08.08.rds')
model_performance = read_rds('./data/model.outputs.2024-08-08.rds')


model_performance[[5]]


str(final_splm)

final_splm$fitted$response
final_splm$obdata$wtmp_mo
library(quantreg)
resid(final_splm)
pred_by_temp = data.frame(WaterTemp = final_splm$obdata$wtmp_mo,
                          COMID = final_splm$obdata$COMID,
                          Year = final_splm$obdata$year,
                          pred = final_splm$fitted$response,
           Resid = abs(resid(final_splm, type = "pearson")),
           Resid_t = resid(final_splm))


sumpred = pred_by_temp %>% 
  summarise(wtmp_mo = mean(WaterTemp), 
            pred = mean(pred),
            Resid_t = mean(Resid_t),
            .by = c(COMID, Year))

multi_rqfit <- rq(wtmp_mo ~ pred, data = sumpred, tau = c(.1,.25,.5,.75,.9))
multi_rqfit
summary(multi_rqfit)
library(emmeans)
test(emtrends(multi_rqfit, ~tau,
                  var = c("pred")), null =1)

library(ggeffects)

quants = data.frame(ggemmeans(multi_rqfit, terms = "pred[all]"))


f1 = ggplot(data = sumpred,
       aes(x = pred,
           y = wtmp_mo)) +
  geom_hex(bins = 100) + 
  geom_ribbon(data = quants, inherit.aes = F,
              aes(x = x, y = predicted,ymin = conf.low, ymax = conf.high,
                  group = group),
              alpha = 0.25)+
  geom_line(data = quants, inherit.aes = F,
            aes(x = x, y = predicted, group = group))+
  annotate("text", label = c("10%", "25%", "50%", "75%", "90%"),
           x = c(7.5,7.5,7,7.25,7.25),
           y = c(3,5.75,8.5,10.5,12), size = 2.5)+
    xlim(2.5, 35.5) + ylim(2.5, 35.5) +
  xlab(expression("Predicted Summer Stream Temperature"~(degree*C))) + 
  ylab(expression("Observed Summer Stream Temperature"~(degree*C))) +
  #ylab('Observed Summer Stream Temperature') +
  scale_fill_gradient(name = 'Count',
                      low = "lightgrey", 
                      high = "red") + 
  geom_abline(color='black', slope = 1, intercept = 0,
              linetype = "dashed") + 
  theme_bw()+
  theme(axis.text = element_text(color = "black"))

library(cowplot)
library(patchwork)


plot_grid(
  plot_grid(maptemps, nrow = 1, ncol = 1, labels = "A"),
  plot_grid(NULL, f1, NULL, nrow = 1, rel_widths = c(0.4, 1, 0.4),
            labels = c("","B",""),
            label_x = 0.05),
  nrow = 2,
  rel_heights = c(1,0.75)
)

ggsave(file = "c:/Users/mmahon/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/Research/USGS Stream Macros/FIshDiversity/Analyses/Fish_Biodiversity_Analyses/Figures/FigureS10.jpg",
       width = 8,
       height = 8.65,
       units = 'in')#,


##test to see whether residuals are biased through time
#Figure indicates minimal directional bias
ggplot(sumpred, aes(x = Year, y = wtmp_mo - pred))+
  geom_violin(aes(group = Year))+
  stat_smooth()

##regression shows non-significant bias in residuals
myr1 = (lm(Resid_t ~ Year, data = sumpred))
summary(myr1)
emmeans::emmeans(myr1, ~Year, at = list(Year = c(1990,1999,2008,2019)))

