setwd('C:/Users/tangw/Desktop/FCS/figure/')
library(ggplot2)
library(reshape2)
library(car)
library(rstatix)
library(tidyverse)
library(ggpubr)
df<-read.csv('C:/Users/tangw/Desktop/FCS/figure/behavior.csv',header=T)
df # 显示结果

data <- gather(df, key = "period", value = "MoCA", -group, -subject)
data2 <- data %>% 
  group_by(group, period) %>% 
  get_summary_stats(MoCA)
#正态性与方差齐性已经验证过

res.aov<-anova_test(data, dv = MoCA, wid = subject, within = period, between = group,effect.size="pes") 
# 组间两两比较
pwc <- data %>%
  group_by(group) %>%
  pairwise_t_test(
    MoCA ~ period, paired = TRUE,
    p.adjust.method = "bonferroni"
  )
pwc <- pwc %>% add_xy_position(x = "period")# 按时间分列
pwc
pwc$y.position[1]<-27.5##修改星号位置
pwc$y.position[2]<-28##修改星号位置
pwc$p.adj.signif[1]<-c("CG: ns")#CG
pwc$p.adj.signif[2]<-c("TG: ***")#TG

pbg <- data %>%
  group_by(period) %>%
  pairwise_t_test(
    MoCA ~ group, paired = F,
    p.adjust.method = "bonferroni"
  )
pbg <- pbg %>% add_xy_position(x = "group")# 按时间分列
pbg

data3<-rbind(data2[2,],data2[1,],data2[4,],data2[3,])
data3$period<-factor(data3$period,levels=c("pre","post"))
fig1<-ggplot(data = data3, aes(x = period, y = mean, color = group, group = group))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width = 0.2,linewidth=1)+
  geom_line(linewidth=1.5)+
  geom_point(data=data3,mapping=aes(group=group, y=mean,fill=group,shape = group),size=3.5)+
  stat_pvalue_manual(pwc, tip.length = 0, hide.ns = F,label.size = 6,bracket.size = 0.6)+
  scale_fill_manual(values = c("TG" = "#ec0000", "CG" = "#00468b"))+
  scale_color_manual(values = c("TG" = "#ec0000", "CG" = "#00468b"))+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(size=18),
        plot.subtitle = element_text(size=14),
        text = element_text(size = 20,color = 'black'),
        axis.text.x = element_text(size = 18,color = 'black'),
        axis.text.y = element_text(size = 18,color = 'black'),
        axis.line.x = element_line(color="black",linewidth = 1),
        axis.line.y = element_line(color="black",linewidth = 1),
        axis.ticks.length = unit(0.3,'cm'),
        axis.ticks.x = element_line(color = "black",linewidth =1),
        axis.ticks.y = element_line(color = "black",linewidth =1),
        legend.key = element_blank())+
  ggtitle("  The effect of intervention
  on cognitive performance")+
  xlab("Period")+
  ylab("scores of MoCA")+
  scale_y_continuous(expand = expansion(0),limits = c(25,28.5),
                     breaks = seq(25,28,0.5))+
  labs(subtitle = get_test_label(res.aov[3,], detailed = TRUE),text = element_text(size = 10,color = 'black'))
  fig1
  ggsave("rmanova.jpeg", fig1 , width = 5, height = 9, dpi = 600)