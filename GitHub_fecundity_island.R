library(tidyverse)
library(lme4)

## Import data
Data<- read.csv("Data_Zanzibar_git.csv")

## Data wrangling
Data<- Data %>%
  mutate(Nfemales= ifelse(is.na(Nfemales), Noadultworms..colony.,Nfemales), 
         age_cat= as.factor(ifelse(Age <= 8,"C1", ifelse(Age >= 9 & Age <= 12, "C2", "C3"))),
         eggs= Intensity..eggs.10ml.,
         Year=as.factor(Year)) %>%
  na.exclude() ## Two individuals are excluded due to missing values in egg counts


#Only data from year 1
Year1<- Data %>%
  filter(Year==1)

S_h<- Year1 %>%
  mutate(n=as.integer(Nfemales),m=Nomiracidia.genotyped)

# Run sibship
source("Run_sibship.R")  

# Mixed effects model   ## Data from Year 1, using Island, age_cat and intervention as fixed effects
re.mod1 <- lmer(log(eggs)~log(expectN)+ Island + age_cat+ Intervention+ (1|Shehia), data=combeggs)

## At baseline, significant difference between Pemba and Unguja
summary(lm(log(eggs)~log(expectN) + Island + age_cat + Intervention, data=combeggs))  


# Fecundity in each island ------------------------------------------------

CIs<-as.data.frame(confint(re.mod1))


CIs[5,1]<- CIs[3,1]+CIs[5,1]
CIs[5,2]<- CIs[3,2]+CIs[5,2]

plotfecY1<- data.frame(mean = c(fixef(re.mod1)[1], fixef(re.mod1)[1]+fixef(re.mod1)[3]),
                       lwrci= c(CIs[3,1], CIs[5,1]),
                       uprci = c(CIs[3,2], CIs[5,2]),
                       Island = c("Pemba", "Unguja"),
                       Year = c("1", "1"))

plotfecY1 <- plotfecY1 %>%
  mutate(mean = exp(mean), lwrci=exp(lwrci), uprci=exp(uprci))


# Year 5 data --------------------------------

Year5<- Data %>%
  filter(Year==5)

S_h<- Year5 %>%
  mutate(n=as.integer(Nfemales),m=Nomiracidia.genotyped)

# Run sibship
source("Run_sibship.R")  

# Mixed effects model
re.mod5 <- lmer(log(eggs)~log(expectN)+ Island + age_cat+ Intervention+ (1|Shehia), data=combeggs)


## At Y5, significant difference between Pemba and Unguja
summary(lm(log(eggs)~log(expectN) + Island + age_cat + Intervention, data=combeggs))  


# Fecundity in each island ------------------------------------------------

CIs<-as.data.frame(confint(re.mod5))


CIs[5,1]<- CIs[3,1]+CIs[5,1]
CIs[5,2]<- CIs[3,2]+CIs[5,2]

plotfecY5<- data.frame(mean = c(fixef(re.mod5)[1], fixef(re.mod5)[1]+fixef(re.mod5)[3]),
                       lwrci= c(CIs[3,1], CIs[5,1]),
                       uprci = c(CIs[3,2], CIs[5,2]),
                       Island = c("Pemba", "Unguja"),
                       Year = c("5", "5"))

plotfecY5 <- plotfecY5 %>%
  mutate(mean = exp(mean), lwrci=exp(lwrci), uprci=exp(uprci))


fec_both<- rbind(plotfecY1, plotfecY5)

