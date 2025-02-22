#####
#figure 2BC
setwd('C:/Users/tangw/Desktop/FCS/figure/')
library(tidyverse)
library(reshape2)
library(ggpubr)
df<-read.csv('C:/Users/tangw/Desktop/FCS/figure/FCSanova.csv',header=T)
fig1<-ggplot(df, aes(x = Group, y = FCSwithin, color = Group)) +
  geom_violin(trim=F,width=1,size=1.2) +
  geom_line(aes(group = paired), color = "grey70", size = 1,alpha=0.1) +
  geom_point(data=df,mapping=aes(group=Group, y=FCSwithin,fill=Group),size=2.5) +
  geom_boxplot(width=0.2,fill='white',size=1,colour='grey50')+
  scale_fill_manual(values = c("#ec0000", "#00468b"))+
  scale_color_manual(values = c("TG" = "#ec0000", "CG" = "#00468b"))+
  stat_compare_means(comparisons = list(c("TG", "CG")), method = "wilcox.test",
                     paired = T, size = 6,label = "p.signif",tip.length = 0,label.y = 5,
                     bracket.size = 1) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        text = element_text(size = 20,color = 'black'),
        axis.text.x = element_text(size = 20,color = 'black'),
        axis.text.y = element_text(size = 20,color = 'black'),
        axis.line.x = element_line(color="black",size = 1),
        axis.line.y = element_line(color="black",size = 1),
        axis.ticks.length = unit(0.3,'cm'),
        axis.ticks.x = element_line(color = "black",linewidth =1),
        axis.ticks.y = element_line(color = "black",linewidth =1))+
  xlab("Group")+
  ylab("Change of within-network FCS")+
  scale_y_continuous(expand = expansion(0),limits = c(-4,7),
                     breaks = seq(-4,6,2))
fig1
ggsave("FCSwithin.jpeg", fig1 , width = 6, height = 4.71, dpi = 600)

fig2<-ggplot(df, aes(x = Group, y = FCSbetween, color = Group)) +
  geom_violin(trim=F,width=1,size=1.2) +
  geom_line(aes(group = paired), color = "grey70", size = 1,alpha=0.1) +
  geom_point(data=df,mapping=aes(group=Group, y=FCSwithin,fill=Group),size=2.5) +
  geom_boxplot(width=0.2,fill='white',size=1,colour='grey50')+
  scale_fill_manual(values = c("#ec0000", "#00468b"))+
  scale_color_manual(values = c("TG" = "#ec0000", "CG" = "#00468b"))+
  stat_compare_means(comparisons = list(c("TG", "CG")), method = "t.test",
                     paired = T, size = 6,label = "p.signif",tip.length = 0,label.y = 9,
                     bracket.size = 1) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        text = element_text(size = 20,color = 'black'),
        axis.text.x = element_text(size = 20,color = 'black'),
        axis.text.y = element_text(size = 20,color = 'black'),
        axis.line.x = element_line(color="black",size = 1),
        axis.line.y = element_line(color="black",size = 1),
        axis.ticks.length = unit(0.3,'cm'),
        axis.ticks.x = element_line(color = "black",linewidth =1),
        axis.ticks.y = element_line(color = "black",linewidth =1))+
  xlab("Group")+
  ylab("Change of between-network FCS")+
  scale_y_continuous(expand = expansion(0),limits = c(-8,11),
                     breaks = seq(-8,10,2))
fig2
ggsave("FCSbetween.jpeg", fig2 , width = 6.5, height = 5.08, dpi = 600)
#####
#figure 6C
setwd('C:/Users/tangw/Desktop/FCS/figure/')
library(ggplot2)
df<-read.csv('C:/Users/tangw/Desktop/FCS/figure/correlation_diff.csv',header=T)
fig6<-ggplot(df,aes(x=Group,y=Moca,color=Group))+
  geom_violin(trim=F,width=1,size=1)+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        text = element_text(size = 16,color = 'black'),
        axis.text.x = element_text(size = 16,color = 'black'),
        axis.text.y = element_text(size = 16,color = 'black'),
        axis.line.x = element_line(color="black",size = 1),
        axis.line.y = element_line(color="black",size = 1),
        axis.ticks.length = unit(0.3,'cm'),
        axis.ticks.x = element_line(color = "black",linewidth =1),
        axis.ticks.y = element_line(color = "black",linewidth =1),
        legend.position = "top")+
  xlab("Group")+
  ylab("Change of Moca")+
  scale_y_continuous(expand = expansion(0),limits = c(-6,9),
                     breaks = seq(-6,8,2))+
  geom_jitter(data=df,mapping=aes(group=Group, y=Moca,fill=Group),size=4,alpha=0.2,
              width=0.35, height=0, shape=21,stroke = 0)+
  #geom_boxplot(width=0.2,fill='white',size=1)+
  scale_fill_manual(values = c("#ec0000", "#00468b"))+
  scale_color_manual(values = c("High" = "#ec0000", "Low" = "#00468b"))+
  geom_boxplot(width=0.2,fill='white',size=1,colour='black')
fig6
ggsave("violin.jpeg", fig6 , width = 4, height = 7, dpi = 600)