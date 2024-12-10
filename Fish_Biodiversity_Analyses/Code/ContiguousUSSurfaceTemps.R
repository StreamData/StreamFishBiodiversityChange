
LandTempsNL = read.csv("./Data/ContUSSummerTemps.csv")

n1 = LandTempsNL %>%
  pivot_longer(cols = July:August,
               names_to = "Month",
               values_to = "Temp") %>%
  filter(Year %in% c(1950:2020)) %>%
  group_by(Year) %>%
  summarize(AvgSummerTemp = mean(Temp))

mean((n1 %>% filter(Year %in% 1950:2000))$AvgSummerTemp)
sd((n1 %>% filter(Year %in% 1950:2000))$AvgSummerTemp)
mean((n1 %>% filter(Year %in% 1990:1994))$AvgSummerTemp)
sd((n1 %>% filter(Year %in% 1990:1994))$AvgSummerTemp)

##no difference between the 50 year average and 5 year average, 0.5F cooler
t.test((n1 %>% filter(Year %in% 1950:2000))$AvgSummerTemp,
       (n1 %>% filter(Year %in% 1990:1994))$AvgSummerTemp)


mean((n1 %>% filter(Year %in% 2000:2020))$AvgSummerTemp)

#72.312
#72.9240 
data.frame(temps = zoo::rollmean(n1$AvgSummerTemp,5),
           Year = 1954:2020) %>%
ggplot(aes(x = Year, y = temps))+
  geom_hline(yintercept = 72.84559, linetype = "dashed")+
  geom_vline(xintercept = 1994, color = "black")+
  scale_x_continuous(breaks = seq(1960,2020,20))+
  ylab("Average summer surface\ntemperature for contiguous US")+
  geom_line() +
  geom_point()+
  theme_bw()+
  theme(axis.text = element_text(color = "black"))
