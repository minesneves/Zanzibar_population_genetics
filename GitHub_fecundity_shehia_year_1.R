library(tidyverse)
library(lme4)
library(data.table)

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


## Run sibship reconstruction
source("Run_sibship.R")

## Mixed effects model
re.mod <- lmer(log(eggs)~log(expectN)+ Island + age_cat+ Intervention+ (1|Shehia), data=combeggs)

# Obtaining mean fecundity and mean intensity of infection per Sheiha --------
beta <-coef(re.mod)
fec_shehias<- setDT(beta$Shehia, keep.rownames = TRUE)[]
fec_shehias_Y1<- fec_shehias[order(fec_shehias$rn),]

