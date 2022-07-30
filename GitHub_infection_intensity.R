## Import libraries
library(lme4)
library(tidyverse)

## Import data
Data<- read.csv("Data_Zanzibar_git.csv")

## Data wrangling
Data<- Data %>%
  mutate(Nfemales= ifelse(is.na(Nfemales), Noadultworms..colony.,Nfemales), 
         age_cat= as.factor(ifelse(Age <= 8,"C1", ifelse(Age >= 9 & Age <= 12, "C2", "C3"))),
         eggs= Intensity..eggs.10ml.,
         Year=as.factor(Year)) %>%
  na.exclude() ## Two individuals are excluded due to missing values in egg counts

## Filter Pemba's data by year
PY1<- Data %>%
  filter(Island=="Pemba" & Year=="1") 

PY5<- Data %>%
  filter(Island=="Pemba" & Year=="5") 


# Pemba Year 1 - mixed effects model -----------------------------------------------
re_modPY1 <- lmer(log(eggs+1) ~ age_cat + (1|Shehia), data=PY1)
AIC(re_modPY1)

fixef_PY1<- fixef(re_modPY1)

CIPY1<-as.data.frame(confint(re_modPY1))

CIPY1[4,1]<- CIPY1[3,1]+CIPY1[4,1]
CIPY1[4,2]<- CIPY1[3,2]+CIPY1[4,2]
CIPY1[5,1]<- CIPY1[3,1]+CIPY1[5,1]
CIPY1[5,2]<- CIPY1[3,2]+CIPY1[5,2]

# Pemba Year 5 - mixed effects model -----------------------------------------------
re_modPY5 <- lmer(log(eggs+1) ~ age_cat + (1|Shehia), data=PY5)

fixef_PY5<- fixef(re_modPY5)

CIPY5<-as.data.frame(confint(re_modPY5))

CIPY5[4,1]<- CIPY5[3,1]+CIPY5[4,1]
CIPY5[4,2]<- CIPY5[3,2]+CIPY5[4,2]

plotfecPY5<- data.frame(mean = c(fixef_PY5[1], 
                                 (fixef_PY5[1]+fixef_PY5[2])),
                        lwrci= c(CIPY5[3,1], CIPY5[4,1]),
                        uprci = c(CIPY5[3,2], CIPY5[4,2]),
                        Island = rep("Pemba", 2),
                        Year = rep("5", 2),
                        age_cat= c("C2", "C3"))



# Unguja ------------------------------------------------------------------

UY1<- Data %>%
  filter(Island=="Unguja" & Year=="1") 

UY5<- Data %>%
  filter(Island=="Unguja" & Year=="5") 


# Unguja Year 1 - mixed effects model -----------------------------------------------
re_modUY1 <- lmer(log(eggs+1) ~ age_cat + (1|Shehia), data=UY1)

fixef_UY1<- fixef(re_modUY1)

CIUY1<-as.data.frame(confint(re_modUY1))

CIUY1[4,1]<- CIUY1[3,1]+CIUY1[4,1]
CIUY1[4,2]<- CIUY1[3,2]+CIUY1[4,2]
CIUY1[5,1]<- CIUY1[3,1]+CIUY1[5,1]
CIUY1[5,2]<- CIUY1[3,2]+CIUY1[5,2]

plotfecUY1<- data.frame(mean = c(fixef_UY1[1], 
                                 (fixef_UY1[1]+fixef_UY1[2]),
                                 (fixef_UY1[1]+fixef_UY1[3])),
                        lwrci= c(CIUY1[3,1], CIUY1[4,1], CIUY1[5,1]),
                        uprci = c(CIUY1[3,2], CIUY1[4,2], CIUY1[5,2]),
                        Island = rep("Unguja", 3),
                        Year = rep("1", 3),
                        age_cat = c("C1", "C2", "C3"))

# Unguja Year 5 - mixed effects model -----------------------------------------------
re_modUY5 <- lmer(log(eggs+1) ~ age_cat + (1|Shehia), data=UY5)

fixef_UY5<- fixef(re_modUY5)

CIUY5<-as.data.frame(confint(re_modUY5))

CIUY5[4,1]<- CIUY5[3,1]+CIUY5[4,1]
CIUY5[4,2]<- CIUY5[3,2]+CIUY5[4,2]

plotfecUY5<- data.frame(mean = c(fixef_UY5[1], 
                                 (fixef_UY5[1]+fixef_UY5[2])),
                        lwrci= c(CIUY5[3,1], CIUY5[4,1]),
                        uprci = c(CIUY5[3,2], CIUY5[4,2]),
                        Island = rep("Unguja", 2),
                        Year = rep("5", 2),
                        age_cat= c("C2", "C3"))

# Checking all mixed effects models used ------------------------------------------------

summary(re_modPY1)
anova(re_modPY1)
summary(re_modPY5)
anova(re_modPY5)
summary(re_modUY1)
anova(re_modUY1)
summary(re_modUY5)
anova(re_modUY5)

# Checking for differences in infection intensity by year in both islands -------------------------------------

summary(lm(log(eggs+1) ~ Island, data=Data))  ## Significant difference between islands across all years and age-groups

summary(lm(log(eggs+1) ~ Island + Year, data=Data))  ## No sign diff between Y1 and Y5 on both islands

summary(lm(log(eggs+1) ~ Island + Year + age_cat, data=Data))  ## No sign diff between Y1 and Y5 but significant differences in age categories on both islands

## Logistic regression for prevalence of heavy infections
Data_log_r<- Data %>% mutate(heavy_inf=ifelse(eggs>=50,1,0))

summary(glm(formula = heavy_inf ~ Island, family = binomial(link = "logit"), 
            data = Data_log_r)) ## across all years and age-groups

summary(glm(formula = heavy_inf ~ Island+Year, family = binomial(link = "logit"), 
            data = Data_log_r)) ## across age-groups



# Checking for differences in infection intensity by age category on each island separately -------------------------------------

# Children 6-8 ------------------------------------------------------------

Y1_C1<- Data %>%
  filter(Year=="1" & age_cat=="C1")

Y1_C1$Island<- fct_relevel(Y1_C1$Island, "Unguja", "Pemba")

summary(lm(log(eggs+1) ~ Island, data=Y1_C1)) ## Significant differences between pemba & unguja in Y1 and children 9-12


# Children 9-12 ------------------------------------------------------------


Y1_C2<- Data %>%
  filter(Year=="1" & age_cat=="C2")

Y1_C2$Island<- fct_relevel(Y1_C2$Island, "Unguja", "Pemba")

Y5_C2<- Data %>%
  filter(Year=="5" & age_cat=="C2")

Y5_C2$Island<- fct_relevel(Y5_C2$Island, "Unguja", "Pemba")


summary(lm(log(eggs+1) ~ Island, data=Y1_C2)) ## Significant differences between pemba & unguja in Y1 and children 9-12

summary(lm(log(eggs+1) ~ Island, data=Y5_C2)) ## Significant differences between pemba & unguja in Y5 and children 9-12


# Adults ------------------------------------------------------------------
Y1_C3<- Data %>%
  filter(Year=="1" & age_cat=="C3")

Y1_C3$Island<- fct_relevel(Y1_C3$Island, "Unguja", "Pemba")

Y5_C3<- Data %>%
  filter(Year=="5" & age_cat=="C3")

Y5_C3$Island<- fct_relevel(Y5_C3$Island, "Unguja", "Pemba")



summary(lm(log(eggs+1) ~ Island, data=Y1_C3)) ## No significant differences between pemba & unguja in Y1 and adults

summary(lm(log(eggs+1) ~ Island, data=Y5_C3)) ## No significant differences between pemba & unguja in Y5 and adults


