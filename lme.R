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