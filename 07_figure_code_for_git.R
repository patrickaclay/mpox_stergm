library(here)
library(ggplot2)
library(GGally)
library(cowplot)
library(dplyr)

############
####Difference in cases averted- primary 
first_pri_cum <- read.csv("data_sheets/cum_cases_best_1dpri.csv")
inter_cum <- read.csv("data_sheets/cum_cases_best_int.csv")
second_pri_cum <- read.csv("data_sheets/cum_cases_best_2dpri.csv")


nyc_data <- read.csv("data_sheets/NYC_data_through_may16.csv")
nyc_data$Date.to.use <- as.Date(nyc_data$Date.to.use, "%Y-%m-%d")

#subtract cumulative cases with first priority from cumulative
#cases with second priority
#So if positive, intermediate is worse
#if negative, intermediate is better
inter_cum_diff <- inter_cum - first_pri_cum

inter_cum_diff_summ <- data.frame(inc.med = numeric(length = 362),
                                  inc.IQR1 = numeric(length = 362),
                                  inc.IQR3 = numeric(length = 362))
for(i in 2:362){
  inter_cum_diff_summ$inc.med[i] <- summary(as.numeric(inter_cum_diff[i,]))[3]
  inter_cum_diff_summ$inc.IQR1[i] <- summary(as.numeric(inter_cum_diff[i,]))[2]
  inter_cum_diff_summ$inc.IQR3[i] <- summary(as.numeric(inter_cum_diff[i,]))[5]
}

#subtract cumulative cases with first priority from cumulative
#cases with second priority
#So if positive, 2dpri is worse
#if negative, 2dpri is better
second_pri_cum_diff <- second_pri_cum - first_pri_cum

second_pri_cum_diff_summ <- data.frame(inc.med = numeric(length = 362),
                                       inc.IQR1 = numeric(length = 362),
                                       inc.IQR3 = numeric(length = 362))
for(i in 2:362){
  second_pri_cum_diff_summ$inc.med[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[3]
  second_pri_cum_diff_summ$inc.IQR1[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[2]
  second_pri_cum_diff_summ$inc.IQR3[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[5]
}




inter_cum_diff_summ$Date.to.use <- nyc_data$Date.to.use
second_pri_cum_diff_summ$Date.to.use <- nyc_data$Date.to.use


second_pri_cum_diff_summ$Strategy <- "Second-dose priority"
inter_cum_diff_summ$Strategy <- "Intermediate"

vaccine_impact_diff <- rbind(second_pri_cum_diff_summ,inter_cum_diff_summ)
vaccine_impact_diff$Strategy <- as.factor(vaccine_impact_diff$Strategy)


scenarios_compare_diff <-   ggplot(data = vaccine_impact_diff, aes(color = Strategy, fill = Strategy)) +
  theme_bw() +
  theme(panel.grid.major = element_line(linewidth = 1), panel.grid.minor = element_line(linewidth = 0.5)) + 
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x = element_text(size=12)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y = element_text(size=12)) +  #size of y-axis title
  theme(legend.text=element_text(size=12)) +
  theme(title =element_text(size=16, face='bold'), plot.title = element_text(hjust = 0.5))  +    #Size of legend Title
  geom_ribbon(aes(ymin = inc.IQR1, ymax = inc.IQR3, x = Date.to.use), alpha = 0.3) +
  geom_line(aes(y = inc.med, x = Date.to.use),linewidth = 2) +
  labs(x = "Date",y="Change in Cases Compared \nto First-Dose Priority ") +
  scale_x_date(date_labels = "%b %Y") +
  scale_color_manual(values = c("darkblue","darkred")) + 
  scale_fill_manual(values = c("darkblue","darkred")) + 
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  ylim(-600,3200)
#geom_text(mapping = aes(x = -Inf, y = -Inf), label = doses$Label,hjust = -0.1,vjust = -1)


ggsave("figures/difference_in_averted_primary.png", scenarios_compare_diff, bg = "transparent",
       width = 7.225, height = 4, dpi = 300, units = "in")



############
####Difference in cases averted- limited
first_pri_cum <- read.csv("data_sheets/cum_cases_limited_1dpri.csv")
inter_cum <- read.csv("data_sheets/cum_cases_limited_int.csv")
second_pri_cum <- read.csv("data_sheets/cum_cases_limited_2dpri.csv")


nyc_data <- read.csv("data_sheets/NYC_data_through_may16.csv")
nyc_data$Date.to.use <- as.Date(nyc_data$Date.to.use, "%Y-%m-%d")


inter_cum_diff <- inter_cum - first_pri_cum

inter_cum_diff_summ <- data.frame(inc.med = numeric(length = 362),
                                  inc.IQR1 = numeric(length = 362),
                                  inc.IQR3 = numeric(length = 362))
for(i in 2:362){
  inter_cum_diff_summ$inc.med[i] <- summary(as.numeric(inter_cum_diff[i,]))[3]
  inter_cum_diff_summ$inc.IQR1[i] <- summary(as.numeric(inter_cum_diff[i,]))[2]
  inter_cum_diff_summ$inc.IQR3[i] <- summary(as.numeric(inter_cum_diff[i,]))[5]
}

second_pri_cum_diff <- second_pri_cum - first_pri_cum

second_pri_cum_diff_summ <- data.frame(inc.med = numeric(length = 362),
                                       inc.IQR1 = numeric(length = 362),
                                       inc.IQR3 = numeric(length = 362))
for(i in 2:362){
  second_pri_cum_diff_summ$inc.med[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[3]
  second_pri_cum_diff_summ$inc.IQR1[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[2]
  second_pri_cum_diff_summ$inc.IQR3[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[5]
}




inter_cum_diff_summ$Date.to.use <- nyc_data$Date.to.use
second_pri_cum_diff_summ$Date.to.use <- nyc_data$Date.to.use


second_pri_cum_diff_summ$Strategy <- "Second-dose priority"
inter_cum_diff_summ$Strategy <- "Intermediate"

vaccine_impact_diff <- rbind(second_pri_cum_diff_summ,inter_cum_diff_summ)
vaccine_impact_diff$Strategy <- as.factor(vaccine_impact_diff$Strategy)


scenarios_compare_diff <-   ggplot(data = vaccine_impact_diff, aes(color = Strategy, fill = Strategy)) +
  theme_bw() +
  theme(panel.grid.major = element_line(linewidth = 1), panel.grid.minor = element_line(linewidth = 0.5)) + 
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x = element_text(size=12)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y = element_text(size=12)) +  #size of y-axis title
  theme(legend.text=element_text(size=12)) +
  theme(title =element_text(size=16, face='bold'), plot.title = element_text(hjust = 0.5))  +    #Size of legend Title
  geom_ribbon(aes(ymin = inc.IQR1, ymax = inc.IQR3, x = Date.to.use), alpha = 0.3) +
  geom_line(aes(y = inc.med, x = Date.to.use),linewidth = 2) +
  labs(x = "Date",y="Change in Cases Compared \nto First-Dose Priority ") +
  scale_x_date(date_labels = "%b %Y") +
  scale_color_manual(values = c("darkblue","darkred")) + 
  scale_fill_manual(values = c("darkblue","darkred")) + 
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
ylim(-600,3200)
#geom_text(mapping = aes(x = -Inf, y = -Inf), label = doses$Label,hjust = -0.1,vjust = -1)


ggsave("figures/difference_in_averted_limited.png", scenarios_compare_diff, bg = "transparent",
       width = 7.225, height = 4, dpi = 300, units = "in")



############
####Difference in cases averted- low incremental effectiveness 
first_pri_cum <- read.csv("data_sheets/cum_cases_best_1dpri_low_inc_VE.csv")
inter_cum <- read.csv("data_sheets/cum_cases_best_int_low_inc_VE.csv")
second_pri_cum <- read.csv("data_sheets/cum_cases_best_2dpri_low_inc_VE.csv")


nyc_data <- read.csv("data_sheets/NYC_data_through_may16.csv")
nyc_data$Date.to.use <- as.Date(nyc_data$Date.to.use, "%Y-%m-%d")


inter_cum_diff <- inter_cum - first_pri_cum

inter_cum_diff_summ <- data.frame(inc.med = numeric(length = 362),
                                  inc.IQR1 = numeric(length = 362),
                                  inc.IQR3 = numeric(length = 362))
for(i in 2:362){
  inter_cum_diff_summ$inc.med[i] <- summary(as.numeric(inter_cum_diff[i,]))[3]
  inter_cum_diff_summ$inc.IQR1[i] <- summary(as.numeric(inter_cum_diff[i,]))[2]
  inter_cum_diff_summ$inc.IQR3[i] <- summary(as.numeric(inter_cum_diff[i,]))[5]
}

second_pri_cum_diff <- second_pri_cum - first_pri_cum

second_pri_cum_diff_summ <- data.frame(inc.med = numeric(length = 362),
                                       inc.IQR1 = numeric(length = 362),
                                       inc.IQR3 = numeric(length = 362))
for(i in 2:362){
  second_pri_cum_diff_summ$inc.med[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[3]
  second_pri_cum_diff_summ$inc.IQR1[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[2]
  second_pri_cum_diff_summ$inc.IQR3[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[5]
}




inter_cum_diff_summ$Date.to.use <- nyc_data$Date.to.use
second_pri_cum_diff_summ$Date.to.use <- nyc_data$Date.to.use


second_pri_cum_diff_summ$Strategy <- "Second-dose priority"
inter_cum_diff_summ$Strategy <- "Intermediate"

vaccine_impact_diff <- rbind(second_pri_cum_diff_summ,inter_cum_diff_summ)
vaccine_impact_diff$Strategy <- as.factor(vaccine_impact_diff$Strategy)


scenarios_compare_diff <-   ggplot(data = vaccine_impact_diff, aes(color = Strategy, fill = Strategy)) +
  theme_bw() +
  theme(panel.grid.major = element_line(linewidth = 1), panel.grid.minor = element_line(linewidth = 0.5)) + 
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x = element_text(size=12)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y = element_text(size=12)) +  #size of y-axis title
  theme(legend.text=element_text(size=12)) +
  theme(title =element_text(size=16, face='bold'), plot.title = element_text(hjust = 0.5))  +    #Size of legend Title
  geom_ribbon(aes(ymin = inc.IQR1, ymax = inc.IQR3, x = Date.to.use), alpha = 0.3) +
  geom_line(aes(y = inc.med, x = Date.to.use),linewidth = 2) +
  labs(x = "Date",y="Change in Cases Compared \nto First-Dose Priority ") +
  scale_x_date(date_labels = "%b %Y") +
  scale_color_manual(values = c("darkblue","darkred")) + 
  scale_fill_manual(values = c("darkblue","darkred")) + 
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  ylim(-600,2700)
#geom_text(mapping = aes(x = -Inf, y = -Inf), label = doses$Label,hjust = -0.1,vjust = -1)


ggsave("figures/difference_in_averted_low_inc_VE.png", scenarios_compare_diff, bg = "transparent",
       width = 7.225, height = 4, dpi = 300, units = "in")


############
####Difference in cases averted- high incremental effectiveness 
first_pri_cum <- read.csv("data_sheets/cum_cases_best_1dpri_high_inc_VE.csv")
inter_cum <- read.csv("data_sheets/cum_cases_best_int_high_inc_VE.csv")
second_pri_cum <- read.csv("data_sheets/cum_cases_best_2dpri_high_inc_VE.csv")


nyc_data <- read.csv("data_sheets/NYC_data_through_may16.csv")
nyc_data$Date.to.use <- as.Date(nyc_data$Date.to.use, "%Y-%m-%d")


inter_cum_diff <- inter_cum - first_pri_cum

inter_cum_diff_summ <- data.frame(inc.med = numeric(length = 362),
                                  inc.IQR1 = numeric(length = 362),
                                  inc.IQR3 = numeric(length = 362))
for(i in 2:362){
  inter_cum_diff_summ$inc.med[i] <- summary(as.numeric(inter_cum_diff[i,]))[3]
  inter_cum_diff_summ$inc.IQR1[i] <- summary(as.numeric(inter_cum_diff[i,]))[2]
  inter_cum_diff_summ$inc.IQR3[i] <- summary(as.numeric(inter_cum_diff[i,]))[5]
}

second_pri_cum_diff <- second_pri_cum - first_pri_cum

second_pri_cum_diff_summ <- data.frame(inc.med = numeric(length = 362),
                                       inc.IQR1 = numeric(length = 362),
                                       inc.IQR3 = numeric(length = 362))
for(i in 2:362){
  second_pri_cum_diff_summ$inc.med[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[3]
  second_pri_cum_diff_summ$inc.IQR1[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[2]
  second_pri_cum_diff_summ$inc.IQR3[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[5]
}




inter_cum_diff_summ$Date.to.use <- nyc_data$Date.to.use
second_pri_cum_diff_summ$Date.to.use <- nyc_data$Date.to.use


second_pri_cum_diff_summ$Strategy <- "Second-dose priority"
inter_cum_diff_summ$Strategy <- "Intermediate"

vaccine_impact_diff <- rbind(second_pri_cum_diff_summ,inter_cum_diff_summ)
vaccine_impact_diff$Strategy <- as.factor(vaccine_impact_diff$Strategy)


scenarios_compare_diff <-   ggplot(data = vaccine_impact_diff, aes(color = Strategy, fill = Strategy)) +
  theme_bw() +
  theme(panel.grid.major = element_line(linewidth = 1), panel.grid.minor = element_line(linewidth = 0.5)) + 
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x = element_text(size=12)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y = element_text(size=12)) +  #size of y-axis title
  theme(legend.text=element_text(size=12)) +
  theme(title =element_text(size=16, face='bold'), plot.title = element_text(hjust = 0.5))  +    #Size of legend Title
  geom_ribbon(aes(ymin = inc.IQR1, ymax = inc.IQR3, x = Date.to.use), alpha = 0.3) +
  geom_line(aes(y = inc.med, x = Date.to.use),linewidth = 2) +
  labs(x = "Date",y="Change in Cases Compared \nto First-Dose Priority ") +
  scale_x_date(date_labels = "%b %Y") +
  scale_color_manual(values = c("darkblue","darkred")) + 
  scale_fill_manual(values = c("darkblue","darkred")) + 
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
ylim(-600,2700)
#geom_text(mapping = aes(x = -Inf, y = -Inf), label = doses$Label,hjust = -0.1,vjust = -1)


ggsave("figures/difference_in_averted_high_inc_VE.png", scenarios_compare_diff, bg = "transparent",
       width = 7.225, height = 4, dpi = 300, units = "in")


############
####Difference in cases averted- clade 1 (higher transmission) 
first_pri_cum <- read.csv("data_sheets/cum_cases_clade1_1dpri.csv")
inter_cum <- read.csv("data_sheets/cum_cases_clade1_int.csv")
second_pri_cum <- read.csv("data_sheets/cum_cases_clade1_2dpri.csv")


nyc_data <- read.csv("data_sheets/NYC_data_through_may16.csv")
nyc_data$Date.to.use <- as.Date(nyc_data$Date.to.use, "%Y-%m-%d")

#subtract cumulative cases with first priority from cumulative
#cases with second priority
#So if positive, intermediate is worse
#if negative, intermediate is better
inter_cum_diff <- inter_cum - first_pri_cum

inter_cum_diff_summ <- data.frame(inc.med = numeric(length = 362),
                                  inc.IQR1 = numeric(length = 362),
                                  inc.IQR3 = numeric(length = 362))
for(i in 2:362){
  inter_cum_diff_summ$inc.med[i] <- summary(as.numeric(inter_cum_diff[i,]))[3]
  inter_cum_diff_summ$inc.IQR1[i] <- summary(as.numeric(inter_cum_diff[i,]))[2]
  inter_cum_diff_summ$inc.IQR3[i] <- summary(as.numeric(inter_cum_diff[i,]))[5]
}

#subtract cumulative cases with first priority from cumulative
#cases with second priority
#So if positive, 2dpri is worse
#if negative, 2dpri is better
second_pri_cum_diff <- second_pri_cum - first_pri_cum

second_pri_cum_diff_summ <- data.frame(inc.med = numeric(length = 362),
                                       inc.IQR1 = numeric(length = 362),
                                       inc.IQR3 = numeric(length = 362))
for(i in 2:362){
  second_pri_cum_diff_summ$inc.med[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[3]
  second_pri_cum_diff_summ$inc.IQR1[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[2]
  second_pri_cum_diff_summ$inc.IQR3[i] <- summary(as.numeric(second_pri_cum_diff[i,]))[5]
}




inter_cum_diff_summ$Date.to.use <- nyc_data$Date.to.use
second_pri_cum_diff_summ$Date.to.use <- nyc_data$Date.to.use


second_pri_cum_diff_summ$Strategy <- "Second-dose priority"
inter_cum_diff_summ$Strategy <- "Intermediate"

vaccine_impact_diff <- rbind(second_pri_cum_diff_summ,inter_cum_diff_summ)
vaccine_impact_diff$Strategy <- as.factor(vaccine_impact_diff$Strategy)


scenarios_compare_diff <-   ggplot(data = vaccine_impact_diff, aes(color = Strategy, fill = Strategy)) +
  theme_bw() +
  theme(panel.grid.major = element_line(linewidth = 1), panel.grid.minor = element_line(linewidth = 0.5)) + 
  theme(axis.title.x = element_text(face="bold", size=16),axis.text.x = element_text(size=12)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=16),axis.text.y = element_text(size=12)) +  #size of y-axis title
  theme(legend.text=element_text(size=12)) +
  theme(title =element_text(size=16, face='bold'), plot.title = element_text(hjust = 0.5))  +    #Size of legend Title
  geom_ribbon(aes(ymin = inc.IQR1, ymax = inc.IQR3, x = Date.to.use), alpha = 0.3) +
  geom_line(aes(y = inc.med, x = Date.to.use),linewidth = 2) +
  labs(x = "Date",y="Change in Cases Compared \nto First-Dose Priority ") +
  scale_x_date(date_labels = "%b %Y") +
  scale_color_manual(values = c("darkblue","darkred")) + 
  scale_fill_manual(values = c("darkblue","darkred")) #+ 
  #theme(
  #  panel.background = element_rect(fill = "transparent",colour = NA),
  #  plot.background = element_rect(fill = "transparent",colour = NA)
  #) #+
  #ylim(-600,3200)
#geom_text(mapping = aes(x = -Inf, y = -Inf), label = doses$Label,hjust = -0.1,vjust = -1)


ggsave("figures/difference_in_averted_primary_clade1.png", scenarios_compare_diff, 
       width = 7.225, height = 4, dpi = 300, units = "in")








###### Figure 1
###### Outline of vaccine strategies

doses <- read.csv("data_sheets/strategy_doses.csv")
doses$Week <- as.Date(doses$Week,"%m/%d/%Y")

#intervention_compare$Intervention <- factor(intervention_compare$Intervention, levels = c("No Intervention"))#,"Vaccination","Behavior Change","Vacc. + Behave."))


strategy_doses <- ggplot(data = doses, aes(x = Week,y = Cumulative, color = Strategy, linetype = Strategy)) +
  theme_bw() +
  theme(panel.grid.major = element_line(linewidth = 1), panel.grid.minor = element_line(linewidth = 0.5)) + 
  theme(axis.title.x = element_text(face="bold", size=14),axis.text.x = element_text(size=8)) +  #Size of x-axis title
  theme(axis.title.y = element_text(face="bold", size=14),axis.text.y = element_text(size=10)) +  #size of y-axis title
  theme(title =element_text(size=14, face='bold'), plot.title = element_text(hjust = 0.5))  +    #Size of legend Title
  geom_line(linewidth = 2) +
  labs(x = "Date",y="Cumulative Doses Given") +
  #ylim(c(0,1080)) +
  scale_x_date(date_labels = "%b %Y")+
  facet_wrap(facets = as.factor(doses$Dose), nrow = 3) + 
  scale_linetype_manual(labels = c("Baseline:\nfirst-dose priority","Counterfactual:\nintermediate","Counterfactual:\nsecond-dose priority"),values = c("solid","dotted","dashed")) +
  scale_color_manual(labels = c("Baseline:\nfirst-dose priority","Counterfactual:\nintermediate","Counterfactual:\nsecond-dose priority"),values = c("darkblue","purple","darkred")) + 
  theme(
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  guides(color = guide_legend(byrow = TRUE), linetype = guide_legend(byrow = TRUE)) +
  theme(legend.spacing.y = unit(0.5, "cm"))
  #geom_text(mapping = aes(x = -Inf, y = -Inf), label = doses$Label,hjust = -0.1,vjust = -1)

ggsave("figures/strategy_doses.png", strategy_doses, bg = "transparent",
       width = 5.5, height = 7, dpi = 300, units = "in")

###### Figure 2
###### model fit and total cases averted


first_pri_inc <- read.csv("data_sheets/inc_cases_best_1dpri.csv")
first_pri_inc_summ <- data.frame(inc.med = numeric(length = 362),
                             inc.IQR1 = numeric(length = 362),
                             inc.IQR3 = numeric(length = 362))
for(i in 2:362){
  first_pri_inc_summ$inc.med[i] <- summary(as.numeric(first_pri_inc[i,]))[3]
  first_pri_inc_summ$inc.IQR1[i] <- summary(as.numeric(first_pri_inc[i,]))[2]
  first_pri_inc_summ$inc.IQR3[i] <- summary(as.numeric(first_pri_inc[i,]))[5]
}


nyc_data <- read.csv("data_sheets/NYC_data_through_may16.csv")
nyc_data$Date.to.use <- as.Date(nyc_data$Date.to.use, "%Y-%m-%d")

first_pri_inc_summ$Date.to.use <- nyc_data$Date.to.use


  
  first_pri_cum <- read.csv("data_sheets/cum_cases_best_1dpri.csv")
  first_pri_cum_summ <- data.frame(inc.med = numeric(length = 362),
                               inc.IQR1 = numeric(length = 362),
                               inc.IQR3 = numeric(length = 362))
  for(i in 2:362){
    first_pri_cum_summ$inc.med[i] <- summary(as.numeric(first_pri_cum[i,]))[3]
    first_pri_cum_summ$inc.IQR1[i] <- summary(as.numeric(first_pri_cum[i,]))[2]
    first_pri_cum_summ$inc.IQR3[i] <- summary(as.numeric(first_pri_cum[i,]))[5]
  }
  
  novacc_cum <- read.csv("data_sheets/cum_cases_best_novacc.csv")
  novacc_cum_summ <- data.frame(inc.med = numeric(length = 362),
                            inc.IQR1 = numeric(length = 362),
                            inc.IQR3 = numeric(length = 362))
  for(i in 2:362){
    novacc_cum_summ$inc.med[i] <- summary(as.numeric(novacc_cum[i,]))[3]
    novacc_cum_summ$inc.IQR1[i] <- summary(as.numeric(novacc_cum[i,]))[2]
    novacc_cum_summ$inc.IQR3[i] <- summary(as.numeric(novacc_cum[i,]))[5]
  }
  

  first_pri_cum_summ$Date.to.use <- nyc_data$Date.to.use
  novacc_cum_summ$Date.to.use <- nyc_data$Date.to.use
  
  
  first_pri_cum_summ$Strategy <- "First-dose priority"
  novacc_cum_summ$Strategy <- "No Vaccination"
  
  vaccine_impact <- rbind(first_pri_cum_summ,novacc_cum_summ)
  vaccine_impact$Strategy <- as.factor(vaccine_impact$Strategy)
  
  first_pri_inc_summ_shift <- first_pri_inc_summ
 # first_pri_inc_summ_shift$Date.to.use <- first_pri_inc_summ_shift$Date.to.use + 10
#  first_pri_inc_summ_shift[,1:3] <- first_pri_inc_summ_shift[,1:3] * 1.6
  
model_fit <- ggplot(data = first_pri_inc_summ_shift[1:250,]) +
    theme_bw() +
    theme(panel.grid.major = element_line(linewidth = 1), panel.grid.minor = element_line(linewidth = 0.5)) + 
    theme(axis.title.x = element_text(face="bold", size=16),axis.text.x = element_text(size=12)) +  #Size of x-axis title
    theme(axis.title.y = element_text(face="bold", size=16),axis.text.y = element_text(size=12)) +  #size of y-axis title
    theme(title =element_text(size=16, face='bold'), plot.title = element_text(hjust = 0.5))  +    #Size of legend Title
    geom_ribbon(aes(ymin = inc.IQR1, ymax = inc.IQR3, x = Date.to.use),color = "darkblue", fill = "darkblue",alpha = 0.3) +
    geom_line(aes(y = inc.med, x = Date.to.use),linewidth = 2, color = "darkblue") +
    geom_point(data = nyc_data[1:250,], aes(x = Date.to.use, y = mean.cases), size = 2) +
    labs(x = "Date",y="Incident Cases") +
    scale_x_date(date_labels = "%b %Y") +
    #scale_color_manual(values = c("darkblue","darkred")) + 
    theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA)
    ) +
  ylim(0,80)
  #geom_text(mapping = aes(x = -Inf, y = -Inf), label = doses$Label,hjust = -0.1,vjust = -1)
  
ggsave("figures/Model_fit.png", model_fit, bg = "transparent",
       width = 6, height = 5, dpi = 300, units = "in")
  
total_averted <- ggplot(data = vaccine_impact, aes(color = Strategy, fill = Strategy)) +
    theme_bw() +
    theme(panel.grid.major = element_line(linewidth = 1), panel.grid.minor = element_line(linewidth = 0.5)) + 
    theme(axis.title.x = element_text(face="bold", size=16),axis.text.x = element_text(size=12)) +  #Size of x-axis title
    theme(axis.title.y = element_text(face="bold", size=16),axis.text.y = element_text(size=12)) +  #size of y-axis title
    theme(legend.text=element_text(size=12)) +
    theme(title =element_text(size=16, face='bold'), plot.title = element_text(hjust = 0.5))  +    #Size of legend Title
    geom_ribbon(aes(ymin = inc.IQR1, ymax = inc.IQR3, x = Date.to.use), alpha = 0.3) +
    geom_line(aes(y = inc.med, x = Date.to.use),linewidth = 2) +
    labs(x = "Date",y="Cumulative Cases") +
    scale_color_manual(values = c("darkblue","black")) + 
    scale_fill_manual(values = c("darkblue","black")) + 
    theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA)
    ) #+
  #geom_text(mapping = aes(x = -Inf, y = -Inf), label = doses$Label,hjust = -0.1,vjust = -1)
  
ggsave("figures/total_averted.png", total_averted, bg = "transparent",
       width = 7.2, height = 5, dpi = 300, units = "in")

####Calculate percent averted
summary(as.numeric((novacc_cum[length(novacc_cum$X),3:102])))

summary(as.numeric((novacc_cum[length(novacc_cum$X),3:102] - first_pri_cum[length(first_pri_cum$X),3:102])/
  novacc_cum[length(novacc_cum$X),3:102]))














