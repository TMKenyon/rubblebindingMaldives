# Kenyon et al. 2024 Ecosphere - Experimental Data


# Libraries ----

library(ggplot2)
library(tidyverse)
library(tidylog)
library(gridExtra)
library(vegan)
library(GGally)
library(ggfortify)
library(emmeans)
library(MASS)
library(ggrepel)
library(mgcv)      #for GAMs
library(gratia)    #for GAM plots
library(broom)     #for tidy output
library(MuMIn)     #for model selection and AICc
library(lubridate) #for processing dates
library(gamm4)
library(splines)
library(glmmTMB)
library(DHARMa) # diagnostics for glmmTMB
library(flextable)



# Data import ----

d1 <- read.csv("Rubble_Binding_Ecology-Experimental_Data.csv", header = T) # csv previously named 'rubble_piles_bind_counts_2.csv'
d1 <- d1 %>% filter(ruined == "no") %>% droplevels() # removing data for grid that was thrown onto the reef crest (W site 3 2-3m Set 1 in month 17)



# Data formatting ----

d1$aspect<-factor(d1$aspect)
d1$site<-factor(d1$site)
d1$depth<-factor(d1$depth)
d1$month<-factor(d1$month)
d1$side<-factor(d1$side)

# Combine depth and aspect

d1$aspect_depth <- paste(d1$aspect, d1$depth, sep="_")

d1$aspect_depth <- factor(d1$aspect_depth,
                          levels = c('Lagoon_2-3m', 'SE_2-3m',
                                     'SE_6-7m', 'W_2-3m', 'W_6-7m'), 
                          labels = c('Shelt_Reef-Flat', 'Shelt_Shal', 
                                     'Shelt_Deep','Exp_Shal', 'Exp_Deep'))
levels(d1$aspect_depth)
summary(d1$macroalgae_total_incl_peyson)
summary(d1$macroalgae_total_excl_peyson)
summary(d1$peysonnalia)

# Creating a dataset where the number of binds per grid per binder is standardised to the number of rubble pieces in the grid
# to get number of binds per binder type, per piece:

d1  <- d1 %>% mutate(cyano_mat_pp = ifelse(cyano_mat==0,0,cyano_mat/rubble_count_start),            
                     turf_algae_pp = ifelse(turf_algae==0,0,turf_algae /rubble_count_start),       
                     cca_minus_HBU_pp = ifelse(cca_minus_HBU==0,0,cca_minus_HBU/rubble_count_start),            
                     macroalgae_total_excl_peyson_pp = ifelse(macroalgae_total_excl_peyson==0,0,macroalgae_total_excl_peyson/rubble_count_start), 
                     macroalgae_total_incl_peyson_pp = ifelse(macroalgae_total_incl_peyson==0,0,macroalgae_total_incl_peyson/rubble_count_start), 
                     peysonnalia_pp = ifelse(peysonnalia==0,0,peysonnalia /rubble_count_start),  
                     sponge_encrust_pp =ifelse(sponge_encrust==0,0,sponge_encrust/rubble_count_start),              
                     ascidian_comb_pp =  ifelse(ascidian_comb==0,0,ascidian_comb/rubble_count_start),              
                     bryozoan_comb_pp = ifelse(bryozoan_comb==0,0,bryozoan_comb/rubble_count_start),               
                     anemone_pp  = ifelse(anemone==0,0,anemone/rubble_count_start),                   
                     coral_pp = ifelse(coral==0,0,coral/rubble_count_start),                       
                     other_incl_HBU_excl_serp_pp = ifelse(other_incl_HBU_excl_serp==0,0,other_incl_HBU_excl_serp/rubble_count_start),   
                     serpulid_pp = ifelse(serpulid==0,0,serpulid/rubble_count_start),
                     other_incl_HBU_serpulid_pp = ifelse(other_incl_HBU_serpulid==0,0,other_incl_HBU_serpulid/rubble_count_start)) # combining other and serpulid due to low numbers

d1_long <- gather(d1, binder, bindspp, cyano_mat_pp,             
                  turf_algae_pp, 
                  cca_minus_HBU_pp,            
                  macroalgae_total_excl_peyson_pp, 
                  macroalgae_total_incl_peyson_pp,
                  peysonnalia_pp,         
                  sponge_encrust_pp,               
                  ascidian_comb_pp,          
                  bryozoan_comb_pp,               
                  anemone_pp,                   
                  coral_pp,                        
                  other_incl_HBU_excl_serp_pp,      
                  serpulid_pp, 
                  other_incl_HBU_serpulid_pp, factor_key=TRUE)

d1_long2 <- gather(d1, binder, totalbindspb, cyano_mat,             
                   turf_algae, 
                   cca_minus_HBU,            
                   macroalgae_total_excl_peyson, 
                   macroalgae_total_incl_peyson,
                   peysonnalia,         
                   sponge_encrust,               
                   ascidian_comb,          
                   bryozoan_comb,               
                   anemone,                   
                   coral,                        
                   other_incl_HBU_excl_serp,      
                   serpulid, 
                   other_incl_HBU_serpulid, factor_key=TRUE)

View(d1_long)
# data is showing the number of binds, of each binder type, per rubble piece for each grid topside and each grid underside at each timepoint.

View(d1_long2)
# data is showing the number of binds, of each binder type, per rubble grid topside and per grid underside at each timepoint.



# Overall binding frequency ----

# Random effect of site should technically be included (1 rep per grid = 2 reps per site, per side, per month)
# However, only 2 per site is low (minimum should be 3) so random effects have been removed to run the models effectively.

# Months 6-18 TOPSIDE & UNDERSIDE ----

# except 0, 2, 2.5 (because only topside was recorded for those, and a lot of zeroes)

d1_3_18 <- d1 %>% filter(month!= "2.5") %>% droplevels()
d1_3_18 <- d1_3_18 %>% filter(month!= "2") %>% droplevels()
d1_3_18 <- d1_3_18 %>% filter(month!= "0") %>% droplevels()
levels(d1_3_18$month)

hist(d1_3_18$total_binds_per_piece) 
hist(log(d1_3_18$total_binds_per_piece+1)) # looks better

levels(d1_3_18$grid_unique) # 30 unique grids

mo3_183 <- glmmTMB(log(total_binds_per_piece+1) ~ aspect_depth +
                     month + side +
                     aspect_depth:month +
                     aspect_depth:side +
                     month:side +
                     (1|site_unique/grid_unique), data=d1_3_18) # grid_ID within site included as random effects.

View(d1_3_18)
car::Anova(mo3_183) # sig. effect of aspectdepth:month, remove non-significant terms

summary(mo3_183)

mo3_183_2 <- update(mo3_183, ~ . - aspect_depth:side)

car::Anova(mo3_183_2)

mo3_183_3 <- update(mo3_183_2, ~ . - month:side)

car::Anova(mo3_183_3)

r.squaredGLMM(mo3_183_3) # R squared is 59%
plot(simulateResiduals(fittedModel = mo3_183_3)) # looks good

write.csv(car::Anova(mo3_183_3), file = "car::Anova(mo3_183_3).csv")

emmeans(mo3_183_3, pairwise ~ aspect_depth | month, type = "response")

pred3_18 <- emmeans(mo3_183_3, pairwise ~ month | aspect_depth, type = "response")

pred3_18_1 <- pred3_18$emmeans %>% as.data.frame()
pred3_18_1.t <- flextable(pred3_18_1) %>% 
  colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
  save_as_docx(path = "pred3_18_1.t.docx")

pred3_18_2 <- pred3_18$contrasts %>% as.data.frame()
pred3_18_2.t <- flextable(pred3_18_2) %>% 
  colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
  save_as_docx(path = "pred3_18_2.t.docx")

pred3_18M <- emmeans(mo3_183_3, pairwise ~ aspect_depth | month, type = "response")

pred3_18_1M <- pred3_18M$emmeans %>% as.data.frame()
pred3_18_1.tM <- flextable(pred3_18_1M) %>% 
  colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
  save_as_docx(path = "pred3_18_1.tM.docx")

pred3_18_2M <- pred3_18M$contrasts %>% as.data.frame()
pred3_18_2.tM <- flextable(pred3_18_2M) %>% 
  colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
  save_as_docx(path = "pred3_18_2.tM.docx")


# Predictions per rubble piece ----

# (not per side per rubble piece as above)

# Have to add up topside + underside binds per grid to get total per piece, rather than total per side per piece.
# This is referenced in the Experimental: Stabilised rubble arrays, Bind density section of the results


d1total <- d1 %>% group_by(aspect, depth, site_unique, grid, grid_unique, month, aspect_depth) %>%
  summarise(total_binds_per_piece_both_side = sum(total_binds_per_piece, na.rm = TRUE))

hist(log(d1total$total_binds_per_piece_both_side+1))

mo_d1total <- glmmTMB(log(total_binds_per_piece_both_side+1) ~ aspect_depth * month +
                   (1|site_unique/grid_unique), data=d1total)

plot(simulateResiduals(fittedModel = mo_d1total))

car::Anova(mo_d1total) # sig. effect of aspectdepth:month

emmeans(mo_d1total, pairwise ~ aspect_depth | month, type = "response") 
# number of binds is double that of models above where side is included.


# Months 3-18 TOPSIDE ONLY ----

# Model for comparing just topside incl. Month 3 (is there difference between Month 3 and 6)

d1_3_18T <- d1 %>% filter(month!= "2") %>% droplevels()
d1_3_18T2 <- d1_3_18T %>% filter(month!= "0") %>% droplevels()

levels(d1_3_18T2$month)

# Now filter for topside

d1_3_18T2 <- d1_3_18T2 %>% filter(side!= "Underside") %>% droplevels()

hist(d1_3_18T2$total_binds_per_piece) 
hist(log(d1_3_18T2$total_binds_per_piece+1)) # looks better


mo3_18T2 <- glmmTMB(log(total_binds_per_piece+1) ~ aspect_depth * month +
                      (1|site_unique/grid_unique), data=d1_3_18T2)

plot(simulateResiduals(fittedModel = mo3_18T2)) # looks good

car::Anova(mo3_18T2) # sig. effect of aspectdepth:month

r.squaredGLMM(mo3_18T2) # R squared is 80%

write.csv(car::Anova(mo3_18T2), file = "car::Anova(mo3_18T2).csv")


predmo3_18T2 <- emmeans(mo3_18T2, pairwise ~ month | aspect_depth, type = "response")

predmo3_18T2_1 <- predmo3_18T2$emmeans %>% as.data.frame()
predmo3_18T2_1.t <- flextable(predmo3_18T2_1) %>% 
  colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
  save_as_docx(path = "predmo3_18T2_1.t.docx")

predmo3_18T2_2 <- predmo3_18T2$contrasts %>% as.data.frame()
predmo3_18T2_2.t <- flextable(predmo3_18T2_2) %>% 
  colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
  save_as_docx(path = "predmo3_18T2_2.t.docx")

predmo3_18T3 <- emmeans(mo3_18T2, pairwise ~ aspect_depth | month, type = "response")

predmo3_18T3_1 <- predmo3_18T3$emmeans %>% as.data.frame()
predmo3_18T3_1.t <- flextable(predmo3_18T3_1) %>% 
  colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
  save_as_docx(path = "predmo3_18T3_1.t.docx")

predmo3_18T3_2 <- predmo3_18T3$contrasts %>% as.data.frame()
predmo3_18T3_2.t <- flextable(predmo3_18T3_2) %>% 
  colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
  save_as_docx(path = "predmo3_18T3_2.t.docx")


# Figure 5 (total number binds per piece) ----

datsum <- d1 %>%
  dplyr::group_by(month, aspect_depth, side) %>%
  summarise(meanpp = mean(total_binds_per_piece),
            sdpp = sd(total_binds_per_piece),
            sepp = sd(total_binds_per_piece)/sqrt(n()),
            ci95 = 1.96*sepp)

datsumB <- d1 %>%
  dplyr::group_by(month, aspect_depth) %>%
  summarise(meanpp = mean(total_binds_per_piece),
            sdpp = sd(total_binds_per_piece),
            sepp = sd(total_binds_per_piece)/sqrt(n()),
            ci95 = 1.96*sepp)

View(datsumB) # not working for months 2 and 2.5

levels(d1$aspect_depth)

colsT <- c("Shelt_Reef-Flat" = "#999999", 
           "Shelt_Shal" = "#56B4E9", 
           "Shelt_Deep" = "#0072B2", 
           "Exp_Shal" = "#E69F00", 
           "Exp_Deep" = "#D55E00")

shapesT <- c("Shelt_Reef-Flat" = 19, 
             "Shelt_Shal" = 19, 
             "Shelt_Deep" = 25, 
             "Exp_Shal" = 19, 
             "Exp_Deep" = 25)

legend_titleT <- "Exposure & Depth"

pd <- position_dodge(0.5)

bindsppG <- ggplot(datsum, aes(y=meanpp, x=month, colour = aspect_depth, shape = aspect_depth, fill = aspect_depth)) +
  geom_point(position=pd, size =3) + 
  facet_wrap(~side) + 
  #geom_segment(aes(colour=aspect_depth)) + # not working...
  geom_errorbar(aes(ymin = meanpp-ci95, ymax = meanpp+ci95), width=.2, position = pd) +
  geom_line(aes(x = month, y = meanpp, colour = aspect_depth, 
                group = aspect_depth), linetype = 2, position = pd) +
  labs(x="Time (months)", y = "Number of binds per rubble piece") +
  scale_y_continuous(limits=c(0,3), breaks = c(0,0.5,1.0,1.5,2.0,2.5,3)) +
  scale_colour_manual(legend_titleT, values=colsT, labels=c("Sheltered Reef Flat", "Sheltered Shallow", "Sheltered Deep","Exposed Shallow", "Exposed Deep")) +
  scale_fill_manual(legend_titleT, values = colsT, labels=c("Sheltered Reef Flat", "Sheltered Shallow", "Sheltered Deep","Exposed Shallow", "Exposed Deep")) +
  scale_shape_manual(legend_titleT, values = shapesT, labels=c("Sheltered Reef Flat", "Sheltered Shallow", "Sheltered Deep","Exposed Shallow", "Exposed Deep")) +
  theme_classic() +
  theme(text = element_text(size=13, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position=c("right")) #0.75,0.74

bindsppG


# Binder community composition ----

# Figure S3 & graphical abstract

levels(d1_long$binder)

d1_long$binder <- factor(d1_long$binder,
                      levels = c("cyano_mat_pp",
                                 "turf_algae_pp",
                                 "macroalgae_total_excl_peyson_pp",
                                 "cca_minus_HBU_pp", 
                                 "peysonnalia_pp",               
                                 "sponge_encrust_pp",             
                                 "ascidian_comb_pp",
                                 "bryozoan_comb_pp",
                                 "anemone_pp",
                                 "coral_pp", 
                                 "serpulid_pp", 
                                 "other_incl_HBU_excl_serp_pp"),    
                      labels = c("Cyanobacteria",
                                 "Turf algae",
                                 "Macroalgae",
                                 "Coralline algae", 
                                 "Peyssonnelia",               
                                 "Sponge",             
                                 "Ascidian",
                                 "Bryozoan",
                                 "Anemone",
                                 "Hard coral", 
                                 "Serpulid worm", 
                                 "Other"))

bindcols <- c("Cyanobacteria" = "#427564",
              "Turf algae" = "darkseagreen1",
              "Macroalgae" = "aquamarine3",
              "Coralline algae" =  "lightpink",
              "Peyssonnelia" = "sandybrown",           
              "Sponge" = "khaki1",     
              "Ascidian" = "mediumpurple1",
              "Bryozoan" =   "lightblue2",
              "Anemone" = "indianred1",
              "Coral" = "peachpuff3",
              "Serpulid" = "ivory",
              "Other" =  "ivory3")

pastelcols <- c("Cyanobacteria" = "#427564",
                "Turf algae" = "#E0EEE0",
                "Macroalgae" = "#9BCD9B",
                "Coralline algae" =  "#EEB4B4",
                "Peyssonnelia" = "peachpuff3",           
                "Sponge" = "#F0E68C",     
                "Ascidian" = "#D8BFD8",
                "Bryozoan" =   "lightblue2",
                "Anemone" = "#FA8072",
                "Hard coral" = "peachpuff3",
                "Serpulid worm" = "#EEE8CD",
                "Other" =  "#C7C7C7")

pd <- position_dodge(0.5)

# Month 2.5

countdat2.5 <- d1_long %>% filter(month == "2.5") %>% droplevels() %>%
  group_by(month, aspect_depth, side, binder) %>%
  summarise(mean = mean(bindspp),
            sd = sd(bindspp),
            se = sd/(sqrt(length(bindspp))),
            ci95 = 1.96*se) %>% as.data.frame()

head(countdat2.5)

(countdatplot2.5 <- ggplot(countdat2.5, aes(x=aspect_depth, 
                                           y= mean,
                                           fill=binder)) +
  geom_bar(position="stack", stat = "identity", colour = "black") +  # for Fig. S3, change position = "stack" to position = "fill"
  #geom_text(aes(label=mean), position=position_dodge(width=0.9), vjust=-0.25) +
  labs(x="Habitat", y = "Mean number binds per piece") + 
  scale_fill_manual(name = "Binding type", values= bindcols) +
  facet_wrap(~side) +
  scale_x_discrete(labels = c("Reef flat 2-3m", "Sheltered 2-3m", "Sheltered 6-7m",
                              "Exposed 2-3m", "Exposed 6-7m")) +
  theme_classic() +
  theme(text = element_text(size=20, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "right") +
  theme(axis.text.x.bottom = element_text(size=20,vjust = 0.6,hjust = 0.5,angle = 45)) +
  ggtitle("Month 2.5")) # remove underside as underside was not monitored in month 2.5

# Month 6

countdat6 <- d1_long %>% filter(month == "6") %>% droplevels() %>%
  group_by(month, aspect_depth, side, binder) %>%
  summarise(mean = mean(bindspp),
            sd = sd(bindspp),
            se = sd/(sqrt(length(bindspp))),
            ci95 = 1.96*se) %>% as.data.frame()

(countdatplot6 <- ggplot(countdat6, aes(x=aspect_depth, 
                                       y= mean,
                                       fill=binder)) +
  geom_bar(position="stack", stat = "identity", colour = "black") +  # for Fig. S3, change position = "stack" to position = "fill"
  #geom_text(aes(label=mean), position=position_dodge(width=0.9), vjust=-0.25) +
  labs(x="Habitat", y = "Mean number binds per piece") + 
  scale_fill_manual(name = "Binding organism", values= bindcols) +
  facet_wrap(~side) +
  scale_x_discrete(labels = c("Lag_Shal", "Shelt_Shal", "Shelt_Deep",
                              "Exp_Shal", "Exp_Deep")) +
  theme_classic() +
  theme(text = element_text(size=13, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "none") +
  theme(axis.text.x.bottom = element_text(size=12,vjust = 0.6,hjust = 0.5,angle = 45)) +
  ggtitle("Month 6"))


# Month 10

countdat10 <- d1_long %>% filter(month == "10") %>% droplevels() %>%
  group_by(month, aspect_depth, side, binder) %>%
  summarise(mean = mean(bindspp),
            sd = sd(bindspp),
            se = sd/(sqrt(length(bindspp))),
            ci95 = 1.96*se) %>% as.data.frame()

(countdatplot10 <- ggplot(countdat10, aes(x=aspect_depth, 
                                         y= mean,
                                         fill=binder)) +
  geom_bar(position="stack", stat = "identity", colour = "black") +  # for Fig. S3, change position = "stack" to position = "fill"
  #geom_text(aes(label=mean), position=position_dodge(width=0.9), vjust=-0.25) +
  labs(x="Habitat", y = "Mean number binds per piece") + 
  scale_fill_manual(name = "Binding organism", values= bindcols) +
  facet_wrap(~side) +
  scale_x_discrete(labels = c("Lag_Shal", "Shelt_Shal", "Shelt_Deep",
                              "Exp_Shal", "Exp_Deep")) +
  theme_classic() +
  theme(text = element_text(size=13, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "none") +
  theme(axis.text.x.bottom = element_text(size=12,vjust = 0.6,hjust = 0.5,angle = 45)) +
  ggtitle("Month 10"))


# Month 18

countdat18 <- d1_long %>% filter(month == "17") %>% droplevels() %>%
  group_by(month, aspect_depth, side, binder) %>%
  summarise(mean = mean(bindspp),
            sd = sd(bindspp),
            se = sd/(sqrt(length(bindspp))),
            ci95 = 1.96*se) %>% as.data.frame()

(countdatplot18 <- ggplot(countdat18, aes(x=aspect_depth, 
                                         y= mean,
                                         fill=binder)) +
  geom_bar(position="stack", stat = "identity", colour = "black") +  # for Fig. S3, change position = "stack" to position = "fill"
  #geom_text(aes(label=mean), position=position_dodge(width=0.9), vjust=-0.25) +
  labs(x="Habitat", y = "Mean number binds per piece") + 
  scale_fill_manual(name = "Binding organism", values= bindcols) +
  facet_wrap(~side) +
  scale_x_discrete(labels = c("Lag_Shal", "Shelt_Shal", "Shelt_Deep",
                              "Exp_Shal", "Exp_Deep")) +
  theme_classic() +
  theme(text = element_text(size=13, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "none") +
  theme(axis.text.x.bottom = element_text(size=12,vjust = 0.6,hjust = 0.5,angle = 45)) +
  ggtitle("Month 18"))



# Permanovas ----

# Conducted in PRIMER


topsidedat <- read.csv("topsidedat.csv", header = TRUE)

# Check levels of replication

topsidestat <- topsidedat %>% group_by(aspect_depth, month) %>%
  summarise(cyano_mat_st = n(),
            turf_algae_st = n(),
            cca_minus_HBU_st = n(), 
            macroalgae_total_excl_peyson_st = n(),
            macroalgae_total_incl_peyson_st = n(),
            peysonnalia_st = n(),
            sponge_encrust_st = n(), 
            ascidian_colonial_st = n(),
            ascidian_tunicate_st = n(),
            ascidian_comb_st = n(),
            bryozoan_encrust_st = n(),
            bryozoan_stalked_st = n(),
            bryozoan_comb_st = n(), 
            anemone_st = n(),
            coral_st = n(), 
            other_incl_HBU_excl_serp_st = n(), 
            other_incl_HBU_serpulid_st = n(),
            serpulid_st = n()) %>% as.data.frame()
#View(topsidestat) # 6 reps per aspect per month per group

topsidestat2 <- topsidedat %>% group_by(aspect_depth, month, site, grid_unique) %>%
  summarise(cyano_mat_st = n(),
            turf_algae_st = n(),
            cca_minus_HBU_st = n(), 
            macroalgae_total_excl_peyson_st = n(),
            macroalgae_total_incl_peyson_st = n(),
            peysonnalia_st = n(),
            sponge_encrust_st = n(), 
            ascidian_colonial_st = n(),
            ascidian_tunicate_st = n(),
            ascidian_comb_st = n(),
            bryozoan_encrust_st = n(),
            bryozoan_stalked_st = n(),
            bryozoan_comb_st = n(), 
            anemone_st = n(),
            coral_st = n(), 
            other_incl_HBU_excl_serp_st = n(), 
            other_incl_HBU_serpulid_st = n(),
            serpulid_st = n()) %>% as.data.frame()
#View(topsidestat2) # 1 rep per grid, and 2 per site. Hence random effects of grid and site cannot be included.

View(topsidedat)

# TOPSIDE ----

# All months together ----

topsidecom <- topsidedat %>% #Makes the community matrix (the response variables)
  dplyr::select(cyano_mat_pp, turf_algae_pp, cca_minus_HBU_pp, macroalgae_total_excl_peyson_pp, peysonnalia_pp,
                sponge_encrust_pp, ascidian_comb_pp, bryozoan_comb_pp, anemone_pp,
                coral_pp, other_incl_HBU_excl_serp_pp, serpulid_pp)

topsidemet <- topsidedat %>% #Makes a meta matrix (the explanatory variables)
  dplyr::select(aspect_depth, aspect, depth, month)

tr.topsidecom <- topsidecom ^ 0.25 # transform the community matrix (fourth root, helps to mitigate the influence of outliers, downgrade them, and suitable if many zeroes).

topsidecom.dist <-vegdist(tr.topsidecom, method = "bray")# community distance matrix with bray-curtis distance - euclidean cannot handle a lot of zeroes, bray-curtis is suitable for abundance-based metrics and can handle many zeroes.

# Homogeneity of variance check

b=betadisper(topsidecom.dist,group=topsidemet[,"aspect_depth"]) #checking homogeneity of variance of your factors. Distance matrix is first argument following by factor of your meta matrix (the factors)
anova(b) #no significance, good
d=betadisper(topsidecom.dist,group=topsidemet[,"month"])
anova(d) #significant difference between months

TukeyHSD(d) #but post hoc shows it's only sig. for 2 comparisons

# Permanova is robust to violations of homogeneity of variance IF the design is balanced, which this is.

# Fit the model

topside.adonis <-adonis2(topsidecom.dist ~ aspect_depth * month, 
                        data = topsidedat, method = "bray", permutations = 999) %>% as.data.frame()
topside.adonis
write.csv(topside.adonis, "topside.adonis.csv")

levels(topsidedat$aspect_depth)
levels(topsidedat$month)

# There is a sig. interaction between month and aspect_depth.
# Will split up into months


# Month 3 ----

topsidedat2.5 <- topsidedat %>% filter(month == 2.5) %>% droplevels()


topsidedat2.5com <- topsidedat2.5 %>% #Makes the community matrix (the response variables) for month 2.5 only
  dplyr::select(cyano_mat_pp, turf_algae_pp, cca_minus_HBU_pp, macroalgae_total_excl_peyson_pp, peysonnalia_pp,
                sponge_encrust_pp, ascidian_comb_pp, bryozoan_comb_pp, anemone_pp,
                coral_pp, other_incl_HBU_excl_serp_pp, serpulid_pp)

topsidedat2.5met <- topsidedat2.5 %>% #Makes the explanatory variable matrix
  dplyr::select(aspect_depth)

tr.topsidedat2.5com <- topsidedat2.5com ^ 0.25

topside2.5.dist <-vegdist(tr.topsidedat2.5com, method = "bray")

topside.adonis2.5 <-adonis2(topside2.5.dist ~ aspect_depth,
                           data = topsidedat2.5, method = "bray", permutations = 999)
topside.adonis2.5 

write.csv(topside.adonis2.5 %>% as.data.frame(), 'topside.adonis2.5.csv')

# sig. effect of aspect_depth
# Conduct pairwise comparison

library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

comps <- pairwiseAdonis::pairwise.adonis(tr.topsidedat2.5com[,1:12],topsidedat2.5met$aspect_depth) 
comps
write.csv(comps %>% as.data.frame(), 'compstop2.5.csv')

# Simper

(sim2_5 <- with(topsidedat2.5met, simper(tr.topsidedat2.5com, aspect_depth, permutations = 999)))
summary(sim2_5)

# Month 6 ----

topsidedat6 <- topsidedat %>% filter(month == 6) %>% droplevels()

topsidedat6com <- topsidedat6 %>% #Makes the community matrix (the response variables) for month 2.5 only
  dplyr::select(cyano_mat_pp, turf_algae_pp, cca_minus_HBU_pp, macroalgae_total_excl_peyson_pp, peysonnalia_pp,
                sponge_encrust_pp, ascidian_comb_pp, bryozoan_comb_pp, anemone_pp,
                coral_pp, other_incl_HBU_excl_serp_pp, serpulid_pp)

topsidedat6met <- topsidedat6 %>% #Makes the explanatory variable matrix
  dplyr::select(aspect_depth)

tr.topsidedat6com <- topsidedat6com ^ 0.25

topside6.dist <-vegdist(tr.topsidedat6com, method = "bray")

topside.adonis6 <-adonis2(topside6.dist ~ aspect_depth,
                         data = topsidedat6, method = "bray", permutations = 999)
topside.adonis6

write.csv(topside.adonis6 %>% as.data.frame(), "topside.adonis6.csv")

comps <- pairwiseAdonis::pairwise.adonis(tr.topsidedat6com[,1:12],topsidedat6met$aspect_depth) 

write.csv(comps %>% as.data.frame(), 'compstop6.csv')

# Simper

(sim6 <- with(topsidedat6met, simper(tr.topsidedat6com, aspect_depth, permutations = 99)))
summary(sim6)

# Month 10 ----

topsidedat10 <- topsidedat %>% filter(month == 10) %>% droplevels()

topsidedat10com <- topsidedat10 %>% #Makes the community matrix (the response variables) for month 2.5 only
  dplyr::select(cyano_mat_pp, turf_algae_pp, cca_minus_HBU_pp, macroalgae_total_excl_peyson_pp, peysonnalia_pp,
                sponge_encrust_pp, ascidian_comb_pp, bryozoan_comb_pp, anemone_pp,
                coral_pp, other_incl_HBU_excl_serp_pp, serpulid_pp)

topsidedat10met <- topsidedat10 %>% #Makes the explanatory variable matrix
  dplyr::select(aspect_depth)

tr.topsidedat10com <- topsidedat10com ^ 0.25

topside10.dist <-vegdist(tr.topsidedat10com, method = "bray")

topside.adonis10 <-adonis2(topside10.dist ~ aspect_depth,
                          data = topsidedat10, method = "bray", permutations = 999)
topside.adonis10

write.csv(topside.adonis10 %>% as.data.frame(), "topside.adonis10.csv")

comps <- pairwiseAdonis::pairwise.adonis(tr.topsidedat10com[,1:12],topsidedat10met$aspect_depth) 
comps
write.csv(comps %>% as.data.frame(), 'compstop10.csv')

# Simper

(sim10 <- with(topsidedat10met, simper(tr.topsidedat10com, aspect_depth, permutations = 99)))
summary(sim10)

# Month 18 ----

topsidedat18 <- topsidedat %>% filter(month == 17) %>% droplevels()

topsidedat18com <- topsidedat18 %>% #Makes the community matrix (the response variables) for month 2.5 only
  dplyr::select(cyano_mat_pp, turf_algae_pp, cca_minus_HBU_pp, macroalgae_total_excl_peyson_pp, peysonnalia_pp,
                sponge_encrust_pp, ascidian_comb_pp, bryozoan_comb_pp, anemone_pp,
                coral_pp, other_incl_HBU_excl_serp_pp, serpulid_pp)

topsidedat18met <- topsidedat18 %>% #Makes the explanatory variable matrix
  dplyr::select(aspect_depth)

tr.topsidedat18com <- topsidedat18com ^ 0.25

topside18.dist <-vegdist(tr.topsidedat18com, method = "bray")

topside.adonis18 <-adonis2(topside18.dist ~ aspect_depth,
                          data = topsidedat18, method = "bray", permutations = 999)
topside.adonis18

write.csv(topside.adonis18 %>% as.data.frame(), "topside.adonis18.csv")

comps <- pairwiseAdonis::pairwise.adonis(tr.topsidedat18com[,1:12],topsidedat18met$aspect_depth) 
comps
write.csv(comps %>% as.data.frame(), 'compstop18.csv')

# Simper

(sim18 <- with(topsidedat18met, simper(tr.topsidedat18com, aspect_depth, permutations = 999)))
summary(sim18)


# UNDERSIDE ----

undersidedat <- read.csv("undersidedat.csv", header = TRUE)

undersidedat$aspect_depth <- as.factor(undersidedat$aspect_depth)
undersidedat$aspect <- as.factor(undersidedat$aspect)
undersidedat$depth <- as.factor(undersidedat$depth)
undersidedat$month <- as.factor(undersidedat$month)

# All months together ----

undersidecom <- undersidedat %>% #Makes the community matrix (the response variables)
  dplyr::select(cyano_mat_pp, turf_algae_pp, cca_minus_HBU_pp, macroalgae_total_excl_peyson_pp, peysonnalia_pp,
                sponge_encrust_pp, ascidian_comb_pp, bryozoan_comb_pp, anemone_pp,
                coral_pp, other_incl_HBU_excl_serp_pp, serpulid_pp)

undersidemet <- undersidedat %>% #Makes a meta matrix (the explanatory variables)
  dplyr::select(aspect_depth, aspect, depth, month)

tr.undersidecom <- undersidecom ^ 0.25 # transform the community matrix (fourth root, helps to mitigate the influence of outliers, downgrade them, and suitable if many zeroes).

undersidecom.dist <-vegdist(tr.undersidecom, method = "bray")# community distance matrix with bray-curtis distance - euclidean cannot handle a lot of zeroes, bray-curtis is suitable for abundance-based metrics and can handle many zeroes.

# Homogeneity of variance check

b=betadisper(undersidecom.dist,group=undersidemet[,"aspect_depth"]) 
anova(b) #sig.
d=betadisper(undersidecom.dist,group=undersidemet[,"month"])
anova(d) #no sig.

TukeyHSD(b) #but post hoc shows it's only sig. for 2 comparisons

# Permanova is robust to violations of homogeneity of variance IF the design is balanced, which this is.

# Fit the model

underside.adonis <-adonis2(undersidecom.dist ~ aspect_depth * month, 
                         data = undersidedat, method = "bray", permutations = 999) %>% as.data.frame()
write.csv(underside.adonis, "underside.adonis.csv")

underside.adonis # sig. interaction

write.csv(underside.adonis %>% as.data.frame(), "underside.adonis.csv")


# Month 6 ----

undersidedat6 <- undersidedat %>% filter(month == 6) %>% droplevels()

undersidedat6com <- undersidedat6 %>% #Makes the community matrix (the response variables) for month 2.5 only
  dplyr::select(cyano_mat_pp, turf_algae_pp, cca_minus_HBU_pp, macroalgae_total_excl_peyson_pp, peysonnalia_pp,
                sponge_encrust_pp, ascidian_comb_pp, bryozoan_comb_pp, anemone_pp,
                coral_pp, other_incl_HBU_excl_serp_pp, serpulid_pp)

undersidedat6met <- undersidedat6 %>% #Makes the explanatory variable matrix
  dplyr::select(aspect_depth)

tr.undersidedat6com <- undersidedat6com ^ 0.25

undersidedat6.dist <-vegdist(tr.undersidedat6com, method = "bray")

undersidedat6.adonis <-adonis2(undersidedat6.dist ~ aspect_depth,
                              data = undersidedat6, method = "bray", permutations = 999)
undersidedat6.adonis

write.csv(undersidedat6.adonis %>% as.data.frame(), "undersidedat6.adonis.csv")

compsund6 <- pairwiseAdonis::pairwise.adonis(tr.undersidedat6com[,1:12],undersidedat6met$aspect_depth) 
compsund6
write.csv(compsund6 %>% as.data.frame(), 'compsund6.csv')

# Simper

(simund6 <- with(undersidedat6met, simper(tr.undersidedat6com, aspect_depth, permutations = 99)))
summary(simund6)

# Month 10 ----

undersidedat10 <- undersidedat %>% filter(month == 10) %>% droplevels()

undersidedat10com <- undersidedat10 %>% #Makes the community matrix (the response variables) for month 2.5 only
  dplyr::select(cyano_mat_pp, turf_algae_pp, cca_minus_HBU_pp, macroalgae_total_excl_peyson_pp, peysonnalia_pp,
                sponge_encrust_pp, ascidian_comb_pp, bryozoan_comb_pp, anemone_pp,
                coral_pp, other_incl_HBU_excl_serp_pp, serpulid_pp)

undersidedat10met <- undersidedat10 %>% #Makes the explanatory variable matrix
  dplyr::select(aspect_depth)

tr.undersidedat10com <- undersidedat10com ^ 0.25

undersidedat10.dist <-vegdist(tr.undersidedat10com, method = "bray")

undersidedat10.adonis <-adonis2(undersidedat10.dist ~ aspect_depth,
                               data = undersidedat10, method = "bray", permutations = 999)
undersidedat10.adonis

write.csv(undersidedat10.adonis %>% as.data.frame(), "undersidedat10.adonis.csv")

compsund10 <- pairwiseAdonis::pairwise.adonis(tr.undersidedat10com[,1:12],undersidedat10met$aspect_depth) 
compsund10
write.csv(compsund10 %>% as.data.frame(), 'compsund10.csv')

# Simper

(simund10 <- with(undersidedat10met, simper(tr.undersidedat10com, aspect_depth, permutations = 99)))
summary(simund10)


# Month 18 ----

undersidedat18 <- undersidedat %>% filter(month == 17) %>% droplevels()

undersidedat18com <- undersidedat18 %>% #Makes the community matrix (the response variables) for month 2.5 only
  dplyr::select(cyano_mat_pp, turf_algae_pp, cca_minus_HBU_pp, macroalgae_total_excl_peyson_pp, peysonnalia_pp,
                sponge_encrust_pp, ascidian_comb_pp, bryozoan_comb_pp, anemone_pp,
                coral_pp, other_incl_HBU_excl_serp_pp, serpulid_pp)

undersidedat18met <- undersidedat18 %>% #Makes the explanatory variable matrix
  dplyr::select(aspect_depth)

tr.undersidedat18com <- undersidedat18com ^ 0.25

undersidedat18.dist <-vegdist(tr.undersidedat18com, method = "bray")

undersidedat18.adonis <-adonis2(undersidedat18.dist ~ aspect_depth,
                               data = undersidedat18, method = "bray", permutations = 999)
undersidedat18.adonis

write.csv(undersidedat18.adonis %>% as.data.frame(), "undersidedat18.adonis.csv")

compsund18 <- pairwiseAdonis::pairwise.adonis(tr.undersidedat18com[,1:12],undersidedat18met$aspect_depth) 
compsund18

write.csv(compsund18 %>% as.data.frame(), 'compsund18.csv')

# Simper

(simund18 <- with(undersidedat18met, simper(tr.undersidedat18com, aspect_depth, permutations = 999)))
summary(simund18)


# nMDS Plots ----

# Topside, All months -----

topsidenmds <- metaMDS(topsidecom, distance = "bray", k=2,trymax=200, autotransform=TRUE,expand=FALSE, plot=FALSE)

topsidenmds$stress # The stress is 0.14. A value below 0.2 is acceptable, stress is 1-R2 value.

stressplot(topsidenmds) # Shows that our non-metric fit R2 is 0.98

topsidenmds.sites.scores<- as.data.frame(scores(topsidenmds, display = 'sites'))
topsidenmds.sites.scores<- data.frame(topsidenmds.sites.scores, topsidemet)
topsidenmds.species.scores<- as.data.frame(scores(topsidenmds, display = 'species'))
topsidenmds.species.scores$Binders<- row.names(topsidenmds.species.scores)

head(topsidenmds.sites.scores)
head(topsidenmds.species.scores)

topsidenmds.species.scores$Binders <- factor(topsidenmds.species.scores$Binders,
                                             levels = c("cyano_mat_pp",
                                                        "turf_algae_pp",
                                                        "macroalgae_total_excl_peyson_pp",
                                                        "cca_minus_HBU_pp", 
                                                        "peysonnalia_pp",               
                                                        "sponge_encrust_pp",             
                                                        "ascidian_comb_pp",
                                                        "bryozoan_comb_pp",
                                                        "anemone_pp",
                                                        "coral_pp", 
                                                        "other_incl_HBU_excl_serp_pp",
                                                        "serpulid_pp"),    
                                             labels = c("Cyanobacteria",
                                                        "Turf algae",
                                                        "Macroalgae",
                                                        "CCA", 
                                                        "Peyssonnelia",               
                                                        "Sponge",             
                                                        "Ascidian",
                                                        "Bryozoan",
                                                        "Anemone",
                                                        "Coral", 
                                                        "Serpulid", 
                                                        "Other"))

topsidenmds.sites.scores$aspectdepthmonth <- paste(topsidenmds.sites.scores$aspect_depth, 
                                                   topsidenmds.sites.scores$month, sep="_")

topsidenmds.sites.scores$aspectdepthmonth <- as.factor(topsidenmds.sites.scores$aspectdepthmonth)

topsidenmds.sites.scores$aspectdepthmonth <- 
  factor(topsidenmds.sites.scores$aspectdepthmonth, 
         levels = c("Exp_Deep_10"  ,"Exp_Deep_17" , "Exp_Deep_2.5" , "Exp_Deep_6" , "Exp_Shal_10", "Exp_Shal_17" ,  
                    "Exp_Shal_2.5" ,"Exp_Shal_6", "Shelt_Deep_10" , "Shelt_Deep_17" , "Shelt_Deep_2.5", "Shelt_Deep_6" , 
                    "Shelt_RF_10" , "Shelt_RF_17" , "Shelt_RF_2.5" , "Shelt_RF_6" , "Shelt_Shal_10" , "Shelt_Shal_17" ,
                    "Shelt_Shal_2.5", "Shelt_Shal_6" ),
         labels = c("ED10","ED18","ED2.5","ED6","ES10","ES18",
                    "ES2.5","ES6","SD10","SD18","SD2.5","SD6",
                    "RF10","RF18","RF2.5","RF6","SS10","SS18",
                    "SS2.5","SS6"))

colsT <- c("SS10" = "#56B4E9",
           "SS18" = "#56B4E9",
           "SS2.5" = "#56B4E9",
           "SS6" = "#56B4E9", 
           "SD10" = "#0072B2",
           "SD18" = "#0072B2",
           "SD2.5" = "#0072B2",
           "SD6" = "#0072B2",
           "ES10" = "#E69F00",
           "ES18" = "#E69F00",
           "ES2.5" = "#E69F00",
           "ES6" = "#E69F00",
           "ED10" = "#D55E00",
           "ED18" = "#D55E00",
           "ED2.5" = "#D55E00",
           "ED6" = "#D55E00", 
           "RF10" = "#999999",
           "RF18" = "#999999",
           "RF2.5" = "#999999",
           "RF2.5" = "#999999")

shapesT <- c("SS10" = 19,
             "SS18" = 19,
             "SS2.5" = 19,
             "SS6" = 19, 
             "SD10" = 25,
             "SD18" = 25,
             "SD2.5" = 25,
             "SD6" = 25,
             "ES10" = 19,
             "ES18" = 19,
             "ES2.5" = 19,
             "ES6" = 19,
             "ED10" = 25,
             "ED18" = 25,
             "ED2.5" = 25,
             "ED6" = 25, 
             "RF10" = 19,
             "RF18" = 19,
             "RF2.5" = 19,
             "RF6" = 19)

centroidsA <- aggregate(cbind(NMDS1, NMDS2) ~ aspectdepthmonth, data = topsidenmds.sites.scores, FUN = mean)
centroidsA <- data.frame(centroidsA) # pulling out just the averages which is what I will plot
legend_title <- "Habitat & Month"

levels(topsidenmds.sites.scores$aspectdepthmonth)

# plot for manuscript

gtopside<-ggplot()+
  geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
  geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')+
  #geom_segment(data=topsideshallnmds.species.scores, aes(y=0, x=0, yend=NMDS2, xend=NMDS1), colour = "grey", arrow=arrow(length=unit(0.3, 'lines')))+ #adding the binding organism drivers
  #geom_text(data=topsideshallnmds.species.scores, aes(y=NMDS2, x=NMDS1, label = Binders, hjust="inward",vjust="inward"), show.legend=FALSE, colour = "grey") +
  geom_point(data=centroidsA, aes(y=NMDS2, x=NMDS1, color=aspectdepthmonth, shape = aspectdepthmonth, fill = aspectdepthmonth), size = 3)+
  scale_colour_manual(legend_title, values=colsT) +
  scale_fill_manual(legend_title, values = colsT) +
  scale_shape_manual(legend_title, values = shapesT) +
  geom_text_repel(data=centroidsA, aes(y=NMDS2, x=NMDS1, label = aspectdepthmonth,
                                       hjust="inward",vjust="inward", color=aspectdepthmonth), show.legend=FALSE) +
  theme_classic() +
  ylab("NMDS2") +
  xlab("NMDS1") +
  ggtitle("Topside") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_vline(xintercept=0, linetype="dashed") + 
  geom_hline(yintercept=0, linetype="dashed") +
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "none") 
gtopside

# nMDS segment plot

gtopsideSEG <-ggplot()+
  geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
  geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')+
  geom_segment(data=topsidenmds.species.scores, aes(y=0, x=0, yend=NMDS2, xend=NMDS1), 
               colour = "grey", arrow=arrow(length=unit(0.3, 'lines'))) + 
  geom_text_repel(data=topsidenmds.species.scores, aes(y=NMDS2, x=NMDS1, 
                                                       #hjust="inward",vjust="inward",
                                                       label = Binders), colour = "gray27") + 
  #geom_point(data=topsideshallnmds.sites.scores, aes(y=NMDS2, x=NMDS1), color="white") +
  #geom_point(data=centroidsA, aes(y=NMDS2, x=NMDS1, colour = aspectmonthdepth)) +
  #scale_colour_manual(legend_title, values=cols) +
  theme_classic() +
  ylab("NMDS2") +
  xlab("NMDS1") +
  #ggtitle("Topside") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_vline(xintercept=0, linetype="dashed") + 
  geom_hline(yintercept=0, linetype="dashed") +
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "none") # don't need the legend
gtopsideSEG

grid.arrange (gtopside, gtopsideSEG, nrow = 2)
# 10 x 8.14 inches with and without legend so I can adjust in Illustrator

# PLOT: Underside, All months ----

undersidenmds <- metaMDS(undersidecom, distance = "bray", k=2, trymax=100,
                         autotransform=TRUE, expand=FALSE, plot=FALSE)

undersidenmds$stress # The stress is 0.17. A value below 0.2 is acceptable, stress is 1-R2 value.

stressplot(undersidenmds) # Shows that our non-metric fit R2 is 0.97

undersidenmds.sites.scores<- as.data.frame(scores(undersidenmds, display = 'sites'))
undersidenmds.sites.scores<- data.frame(undersidenmds.sites.scores, undersidemet)
undersidenmds.species.scores<- as.data.frame(scores(undersidenmds, display = 'species'))
undersidenmds.species.scores$Binders<- row.names(undersidenmds.species.scores)

head(undersidenmds.sites.scores)
head(undersidenmds.species.scores)

undersidenmds.species.scores$Binders <- factor(undersidenmds.species.scores$Binders,
                                               levels = c("cyano_mat_pp",
                                                          "turf_algae_pp",
                                                          "macroalgae_total_excl_peyson_pp",
                                                          "cca_minus_HBU_pp", 
                                                          "peysonnalia_pp",               
                                                          "sponge_encrust_pp",             
                                                          "ascidian_comb_pp",
                                                          "bryozoan_comb_pp",
                                                          "anemone_pp",
                                                          "coral_pp", 
                                                          "other_incl_HBU_excl_serp_pp",
                                                          "serpulid_pp"),    
                                               labels = c("Cyanobacteria",
                                                          "Turf algae",
                                                          "Macroalgae",
                                                          "CCA", 
                                                          "Peyssonnelia",               
                                                          "Sponge",             
                                                          "Ascidian",
                                                          "Bryozoan",
                                                          "Anemone",
                                                          "Coral", 
                                                          "Serpulid", 
                                                          "Other"))

undersidenmds.sites.scores$aspectdepthmonth <- paste(undersidenmds.sites.scores$aspect_depth, 
                                                     undersidenmds.sites.scores$month, sep="_")

undersidenmds.sites.scores$aspectdepthmonth <- as.factor(undersidenmds.sites.scores$aspectdepthmonth)

levels(undersidenmds.sites.scores$aspectdepthmonth)

undersidenmds.sites.scores$aspectdepthmonth <- factor(undersidenmds.sites.scores$aspectdepthmonth, 
         levels = c("Exp_Deep_10", "Exp_Deep_17", "Exp_Deep_6", "Exp_Shal_10", "Exp_Shal_17", "Exp_Shal_6", "Shelt_Deep_10",
                    "Shelt_Deep_17", "Shelt_Deep_6" , "Shelt_RF_10", "Shelt_RF_17", "Shelt_RF_6", "Shelt_Shal_10", "Shelt_Shal_17",
                    "Shelt_Shal_6"),
         labels = c("ED10","ED18","ED6","ES10","ES18","ES6","SD10",
                    "SD18","SD6","RF10","RF18","RF6","SS10","SS18",
                    "SS6"))

colsT <- c("SS10" = "#56B4E9",
           "SS18" = "#56B4E9",
           "SS6" = "#56B4E9", 
           "SD10" = "#0072B2",
           "SD18" = "#0072B2",
           "SD6" = "#0072B2",
           "ES10" = "#E69F00",
           "ES18" = "#E69F00",
           "ES6" = "#E69F00",
           "ED10" = "#D55E00",
           "ED18" = "#D55E00",
           "ED6" = "#D55E00", 
           "RF10" = "#999999",
           "RF18" = "#999999",
           "RF6" = "#999999")

shapesT <- c("SS10" = 19,
             "SS18" = 19,
             "SS6" = 19, 
             "SD10" = 25,
             "SD18" = 25,
             "SD6" = 25,
             "ES10" = 19,
             "ES18" = 19,
             "ES6" = 19,
             "ED10" = 25,
             "ED18" = 25,
             "ED6" = 25, 
             "RF10" = 19,
             "RF18" = 19,
             "RF6" = 19)

centroidsA <- aggregate(cbind(NMDS1, NMDS2) ~ aspectdepthmonth, data = undersidenmds.sites.scores, FUN = mean)
centroidsA <- data.frame(centroidsA) # pulling out just the averages which is what I will plot
legend_title <- "Habitat & Month"

# plot for manuscript

gunderside<-ggplot()+
  geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
  geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')+
  #geom_segment(data=topsideshallnmds.species.scores, aes(y=0, x=0, yend=NMDS2, xend=NMDS1), colour = "grey", arrow=arrow(length=unit(0.3, 'lines')))+ #adding the binding organism drivers
  #geom_text(data=topsideshallnmds.species.scores, aes(y=NMDS2, x=NMDS1, label = Binders, hjust="inward",vjust="inward"), show.legend=FALSE, colour = "grey") +
  geom_point(data=centroidsA, aes(y=NMDS2, x=NMDS1, color=aspectdepthmonth, shape = aspectdepthmonth, fill = aspectdepthmonth), size = 3)+
  scale_colour_manual(legend_title, values=colsT) +
  scale_fill_manual(legend_title, values = colsT) +
  scale_shape_manual(legend_title, values = shapesT) +
  geom_text_repel(data=centroidsA, aes(y=NMDS2, x=NMDS1, label = aspectdepthmonth,
                                      hjust="inward",vjust="inward", color=aspectdepthmonth), show.legend=FALSE) +
  theme_classic() +
  ylab("NMDS2") +
  xlab("NMDS1") +
  ggtitle("Underside") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_vline(xintercept=0, linetype="dashed") + 
  geom_hline(yintercept=0, linetype="dashed") +
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "none") 
gunderside

# nMDS segment plot

gundersideSEG <-ggplot()+
  geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
  geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')+
  geom_segment(data=undersidenmds.species.scores, aes(y=0, x=0, yend=NMDS2, xend=NMDS1), 
               colour = "grey", arrow=arrow(length=unit(0.3, 'lines'))) + 
  geom_text_repel(data=undersidenmds.species.scores, aes(y=NMDS2, x=NMDS1, 
                                                         #hjust="inward",vjust="inward",
                                                         label = Binders), colour = "gray27") + 
  #geom_point(data=topsideshallnmds.sites.scores, aes(y=NMDS2, x=NMDS1), color="white") +
  #geom_point(data=centroidsA, aes(y=NMDS2, x=NMDS1, colour = aspectmonthdepth)) +
  #scale_colour_manual(legend_title, values=cols) +
  theme_classic() +
  ylab("NMDS2") +
  xlab("NMDS1") +
  #ggtitle("Topside") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_vline(xintercept=0, linetype="dashed") + 
  geom_hline(yintercept=0, linetype="dashed") +
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "none") # don't need the legend
gundersideSEG

grid.arrange (gunderside, gundersideSEG, nrow = 2)
# 10 x 8.14 inches with and without legend so I can adjust in Illustrator


# Graphical abstract ----

# Prop of rubble bound / Number of binds

# Don't have this in same way as surveys, but assumption is made that each rubble has at least
# 1 bind when the binds pp is 1 (averaged across rubble pieces). This would likely be the case
# as the binds were spread relatiely evenly across the grid.

# Have to add up topside + underside binds per grid to get total per piece, rather than total per side per piece.

d1total <- d1 %>% group_by(aspect, depth, site, grid, grid_unique, month, aspect_depth) %>%
  summarise(total_binds_per_piece_both_side = sum(total_binds_per_piece, na.rm = TRUE))

d1totalsum <- d1total %>%
  dplyr::group_by(month, aspect) %>%
  summarise(meanpp = mean(total_binds_per_piece_both_side),
            sdpp = sd(total_binds_per_piece_both_side),
            sepp = sd(total_binds_per_piece_both_side)/sqrt(n()),
            ci95 = 1.96*sepp)
d1totalsum

# Prop of binds of each binder type (on average)

stats5 <- d1_long %>% 
  dplyr::group_by(month, aspect, side) %>% 
  dplyr::summarise(meanbindspp = sum(bindspp))

stats6 <- d1_long %>% group_by(month, aspect, side, binder) %>% 
  dplyr::summarise(meanbindspp = mean(bindspp)) %>% as.data.frame()
View(stats6)

write.csv(stats6, "stats6.csv")


# Coral cover types before & after bleaching ----


cc1 <- read.csv("Benthic_Data_Maldives_2015-2019.csv", header = T)
cc2 <- read.csv("Benthic_Data_Maldives_2015-2019_2.csv", header = T) 
names(cc2)

# Data formatting ----

cc1$Year<-factor(cc1$Year)
cc1$DepthCat<-factor(cc1$DepthCat)
cc1$TransectID<-factor(cc1$TransectID)
cc1$SiteID<-factor(cc1$SiteID)

cc2$TransectID<-factor(cc2$TransectID)
cc2$CatState<-factor(cc2$CatState)

# Join together ----

cc1_2 <- inner_join(cc1, cc2, by = "TransectID")
str(cc1_2)

# Filter for sites of interest (Vabbin West (8), East (5), South (1)) ----
# North is 2

cc1_2F <- cc1_2 %>% filter(SiteID %in% c("1", "2", "5", "8")) %>% droplevels() # removing data for grid that was thrown onto the reef crest (W site 3 2-3m Set 1 in month 17)
str(cc1_2F)

# Convert to percentages -----

cc1_2F_PC <- cc1_2F %>%
  group_by(TransectID) %>%                       # Group by TransectID
  mutate(Total_Count = sum(Count)) %>%           # Calculate total Count per TransectID
  mutate(Percentage = (Count / Total_Count) * 100) %>%  # Calculate percentage for each CatState
  ungroup()                                      # Ungroup the data

View(cc1_2F_PC)

# Bar plot ----

ggplot(cc1_2F_PC, aes(x = DepthCat, y = Percentage, fill = CatState)) +
  geom_bar(stat = "identity", position = "fill") +  # Bar plot with grouped bars
  facet_grid(SiteID ~ Year) +                              # Separate plots for each Year
  labs(
    title = "Percentage of CatState per Year and Depth",
    x = "CatState",
    y = "Percentage (%)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

# Summary table of benthic composition across sites and years:

cc1_2F_PCsum <- cc1_2F_PC %>% group_by(Year, DepthCat, SiteID, CatState) %>%
  summarise(MeanPercentage = mean(Percentage, na.rm = TRUE),
            sdPercentage = sd(Percentage, na.rm = TRUE),
            sePercentage = sdPercentage/sqrt(n()))
View(cc1_2F_PCsum)
write.csv(cc1_2F_PCsum, "cc1_2F_PCsum.csv")

# % of live coral of different types per site, per depth in each year.

livecoraldat <- cc1_2F_PC %>%
  filter(grepl("L$", CatState)) %>% droplevels()

str(livecoraldat)



ggplot(livecoraldat, aes(x = DepthCat, y = Percentage, fill = CatState)) +
  geom_bar(stat = "identity", position = "stack") +  # Bar plot with grouped bars
  facet_grid(SiteID ~ Year) +                              # Separate plots for each Year
  labs(
    title = "Percentage of CatState per Year and Depth",
    x = "CatState",
    y = "Percentage (%)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels


# Summary table of coral composition across sites and years:

livecoralsum <- livecoraldat %>% group_by(Year, SiteID, DepthCat, CatState) %>%
  summarise(MeanPercentage = mean(Percentage, na.rm = TRUE),
            sdPercentage = sd(Percentage, na.rm = TRUE),
            sePercentage = sdPercentage/sqrt(n()))
View(livecoralsum) # If grouping by SiteID, I looked at Site 8 (West) to write the following passage in the Discussion:
#Reef flat data are not available, but exposed shallow slopes (1-6 m) 
#had 14% cover of digitate corals that fell to 6% post-bleaching; 
#these may have broken down and been more easily transported than 
#branched pieces (Kenyon et al. 2023c) into the reef flat.
#The coordinates for Vabbin West in this data seem to be associated with reef flat,
#but monitoring was for reef slope only... coords wrong.

write.csv(livecoralsum, "livecoralsum.csv")

ggplot(livecoralsum, aes(x = DepthCat, y = MeanPercentage, fill = CatState)) +
  geom_bar(stat = "identity", position = "stack") +  # Bar plot with grouped bars
  facet_wrap(~ Year) +                              # Separate plots for each Year
  labs(
    title = "Percentage of CatState per Year and Depth",
    x = "CatState",
    y = "Percentage (%)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
