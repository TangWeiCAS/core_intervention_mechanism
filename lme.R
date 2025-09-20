setwd('C:/Users/tangw/Desktop/FCS/data/')
library(lme4)
behdata<-read.csv('C:/Users/tangw/Desktop/FCS/correlation_TG2.csv')
#"LLwithin7.csv","RRwithin7.csv","LRwithin7.csv","RLwithin7.csv","LLbetween7.csv","RRbetween7.csv","LRbetween7.csv","RLbetween7.csv"
#"between7","within7"
for(k in c("TG2_within7.csv"))
{FCS<-read.csv(k)
data<-cbind(behdata,FCS)
for(j in c("Moca"))
{
for(i in c("visual","SMN","DA","VA","Limbic","FPCN","DMN"))
  {
  fm <- lm(data[[j]] ~ 1 + age + sex + edu + fd + data[[i]], data)
  print(k)
  print(cbind(j," ~ ",i))
  print(summary(fm))
}
}
}

#####
#全部的被试（加被试随机效应）
setwd('C:/Users/tangw/Desktop/FCS/data/')
library(lme4)
library(lmerTest)
behdata<-read.csv('C:/Users/tangw/Desktop/FCS/correlation.csv')
#"between network FCS"
for(k in c("within network FCS"))
{CG<-read.csv("CG_within7.csv")
TG<-read.csv("TG_within7.csv")
FCS<-rbind(CG,TG)
data<-cbind(behdata,FCS)
for(j in c("Moca"))
{
  for(i in c("visual","SMN","DA","VA","Limbic","FPCN","DMN"))
  {
    fm <- lmer(data[[j]] ~ 1 + age + sex + edu + fd + data[[i]] + (1|subject), data)
    print(k)
    print(cbind(j," ~ ",i))
    print(summary(fm))
  }
}
}

#####
#差值预测改变量
setwd('C:/Users/tangw/Desktop/FCS/data/')
library(lme4)
behdata<-read.csv('C:/Users/tangw/Desktop/FCS/correlation_diff.csv')
for(k in c("diff_within7.csv"))
{FCS<-read.csv(k)
data<-cbind(behdata,FCS)
for(j in c("Moca"))
{
  for(i in c("visual","SMN","DA","VA","Limbic","FPCN","DMN"))
  {
    fm <- lm(data[[j]] ~ 1 + age + sex + edu + data[[i]], data)
    print(k)
    print(cbind(j," ~ ",i))
    print(summary(fm))
  }
}
}

#####
#看看次要指标，说明运动提升身体能力，间接提升情绪水平
#失败了，完全没有，应该是和运动指标直接关联，认知与他无关
setwd('C:/Users/tangw/Desktop/FCS/data/')
library(lme4)
library(lmerTest)
behdata<-read.csv('C:/Users/tangw/Desktop/FCS/correlation.csv')
#"between network FCS"
for(k in c("within network FCS"))
{CG<-read.csv("CG_between7.csv")
TG<-read.csv("TG_between7.csv")
FCS<-rbind(CG,TG)
data<-cbind(behdata,FCS)
for(j in c("SAS","CESD"))
{
  for(i in c("visual","SMN","DA","VA","Limbic","FPCN","DMN"))
  {
    fm <- lmer(data[[j]] ~ 1 + age + sex + edu + fd + data[[i]] + (1|subject), data)
    print(k)
    print(cbind(j," ~ ",i))
    print(summary(fm))
  }
}
}

#####
#所有组差值预测改变量
setwd('C:/Users/tangw/Desktop/FCS/data/')
library(lme4)
library(lmerTest)
data<-read.csv('C:/Users/tangw/Desktop/FCS/correlation_diff2.csv')
for(k in c("diff_within7.csv"))
{
for(j in c("Moca"))
{
  for(i in c("SMN","DA","DMN"))
  {
    fm <- lmer(data[[j]] ~ 1 + age + sex + edu + data[[i]] + (1|group), data)
    print(k)
    print(cbind(j," ~ ",i))
    print(summary(fm))
  }
}
}

#####
#Lme plot (DAN)
library(lme4)
library(ggplot2)
library(ggeffects)
fm <- lmer(Moca ~ 1 + age + sex + edu + fd + DA + (1|subject), data)
summary(fm)
pred.mm<-ggpredict(fm,terms=c("DA [2:25 by=0.1]"))
pred.mm

#绘制预测结果
fig<-ggplot(pred.mm,aes(x, y = predicted)) +
  geom_line(data = data, 
            mapping = aes(x = DA, y = Moca, group = subject, color = group),
            alpha = 1, linewidth = 0.5) +
  geom_point(data,mapping=aes(x = DA, y = Moca, colour = group),size=2,alpha=1) + 
  scale_color_manual(values = c("TG" = "#ec0000", "CG" = "#00468b"))+  
  geom_ribbon(aes(x, ymin = predicted - std.error, ymax = predicted + std.error),      
              fill = "lightgrey", alpha = 0.9) +  
  geom_smooth(method = "lm",alpha=1,linewidth=1.5,colour="black") +
  labs(x = "DAN Within-FCS", y = "MoCA") +  
  theme_bw() +
  theme(panel.grid = element_blank(),legend.position="top",
        panel.border = element_rect(colour = "black", linewidth = 1),
        text = element_text(size = 15),
        axis.title = element_text(size = 14, color = "black",face = "bold"),
        axis.text = element_text(size = 14,color = "black"),
        axis.text.x = element_text(margin=margin(t =3)),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(margin = margin(r = 12)),
        axis.ticks.x = element_line(color = "black",linewidth =0.5),
        axis.ticks.length.x = unit(0.3,"cm"),
        axis.ticks.y = element_line(color = "black",linewidth =0.5),
        axis.ticks.length.y = unit(0.3,"cm"))+
  scale_y_continuous(expand = expansion(0),limits = c(19,31),
                     breaks = seq(20,30,2))+
  scale_x_continuous(expand = expansion(0),limits = c(0,26),
                     breaks = seq(0,26,5))
fig
ggsave("lme.tiff", fig , width = 8, height = 8, dpi = 600)