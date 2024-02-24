# Kenyon et al. (2024) Ecology Submission - Observational Data Code

# Libraries ----

library(glmmTMB)
library(DHARMa)
library(emmeans)
library(tidyverse)
library(MuMIn)
library(Rmisc)
library(flextable)

# Set table output defaults

set_flextable_defaults(
  font.family = "Times New Roman",
  font.size = 10)


# Data ----

dat <- read.csv("/Users/taniakenyon/Dropbox/PhD/Data/Maldives/Binding_Surveys/Binding_Prevalence_Analysis/Rubble_Binding_Ecology-Observational_Data.csv", header=T)
head(dat)
str(dat)
View(dat)



# Formatting ----

cols1 <- c("habitat", "site", "site_unqiue", "depth", "trip", "remove_for_size", "quadrat","square",
           "morphology_simple", "stable", "bound", "binder_1", "binder_2", "binder_3")

dat <- dat %>%  mutate_at(cols1, factor)
levels(dat$bound)

dat <- dat %>% mutate(boundB = ifelse(dat$bound == "Yes", 1, 0)) # creating a column where 'bound' is binary

dat$stable <- factor(dat$stable, levels = c("No","Yes"),
                                 labels = c("Easily mobilised",
                                            "Stable"))
levels(dat$stable)
dat$aspectdepth <- paste(dat$habitat, dat$depth, sep = "_")
dat$aspectdepth <- as.factor(dat$aspectdepth)
dat$widest_span_mm <- as.numeric(dat$widest_span_mm)
class(dat$widest_span_mm)
levels(dat$aspectdepth)



# Data exploration ----

dat2 <- dat %>%
  group_by(aspectdepth, site) %>%
  summarise(count = n_distinct(quadrat)) # 6-9 quadrats per site

ggplot(dat, aes(x =  widest_span_mm, y = boundB, fill = aspectdepth, colour = aspectdepth)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  geom_point(position=position_jitter(width=0.2, height=0)) + 
  facet_wrap(~morphology_simple) +
  theme_classic()

datBR <- dat %>% filter(morphology_simple != "Other")

ggplot(datBR, aes(x =  widest_span_mm/10, y = boundB, fill = morphology_simple, colour = morphology_simple)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  geom_point(position=position_jitter(width=0.2, height=0)) + 
  theme_classic()

ggplot(datBR, aes(x =  widest_span_mm, y = boundB, fill = aspectdepth, colour = aspectdepth)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial")) +
  geom_point(position=position_jitter(width=0.2, height=0)) + 
  theme_classic() # generally, chance of binding is increasing as size increases, in every habitat

# interaction between morphology and size for binding
# if branched, equal chance of being bound across sizes
# if unbranched, low chance of being bound when small, increasing binding chance as size increases

dat.sum2 <- dat %>%
  filter(!is.na(stable)) %>%
  group_by(aspectdepth, morphology_simple, stable) %>%
  dplyr::summarise(meanB = mean(boundB, na.rm=TRUE),
            seB = sd(boundB, na.rm=TRUE)/sqrt(n()),
            N=n())

dat %>% filter(!is.na(bound)) %>% ggplot(aes(aspectdepth)) +
  geom_bar(aes(fill=bound), position = "fill")

# Appears to be less binding in the reef flat, followed by SE, and most binding in the W especially at 6-7m
# Unlikely to be depth effects

dat %>% filter(!is.na(bound)) %>% 
  filter(!is.na(stable)) %>% 
  ggplot(aes(aspectdepth)) +
  geom_bar(aes(fill=bound), position = "fill") +
  facet_grid(stable~morphology_simple)

dat %>% filter(!is.na(bound)) %>% 
  ggplot(aes(morphology_simple)) +
  geom_bar(aes(fill=bound), position = "fill")

# branched pieces bound more often unbranched, particularly when easily mobilised

# Benthic cover (Rubble %) ----

rubcov <- read.csv("~/Dropbox/PhD/Data/Maldives/Benthic_surveys/BenthicSurveys/Benthic_Surveys_FINAL_TM_NAs.csv")

rubcov$habitat <- as.factor(paste(rubcov$Aspect, rubcov$Depth, sep = "_"))

rubcovdat <- summarySE(rubcov, measurevar="percent_rubble", groupvars=c("habitat"), na.rm=TRUE, .drop=FALSE)

rubcovdat # this shows the mean and se for each habitat in terms of rubble cover.

rubcovdat2 <- rubcovdat %>% mutate(cilow = percent_rubble - ci,
                                   cihigh = percent_rubble + ci)

rubcovdat3 <- summarySE(rubcov, measurevar="percent_rubble_less20", groupvars=c("habitat"), na.rm=TRUE, .drop=FALSE)

rubcovdat3 # this shows the mean and se for each habitat in terms of rubble cover.

rubcovdat4 <- rubcovdat3 %>% mutate(cilow = percent_rubble_less20 - ci,
                                   cihigh = percent_rubble_less20 + ci)

rubcovdat5 <- summarySE(rubcov, measurevar="percent_rubble_over20", groupvars=c("habitat"), na.rm=TRUE, .drop=FALSE)

rubcovdat5 # this shows the mean and se for each habitat in terms of rubble cover.

rubcovdat6 <- rubcovdat5 %>% mutate(cilow = percent_rubble_over20 - ci,
                                    cihigh = percent_rubble_over20 + ci)

lm1_1 <- glmmTMB(percent_rubble ~ habitat + (1|Site_Unique/Transect..), 
                 family = gaussian, data = rubcov)

plot(simulateResiduals(fittedModel = lm1_1))

car::Anova(lm1_1)

lm1_1Anova <- car::Anova(lm1_1)

write.csv(car::Anova(lm1_1), "lm1_1Anova.csv")

lm1_1preds <- emmeans(lm1_1, pairwise ~ habitat, type = "response")

write.csv(lm1_1preds$emmeans, "lm1_1predsemmeans.csv")
write.csv(lm1_1preds$contrasts, "lm1_1predscontrasts.csv")

# Rubble bed depth ----

beddat <- summarySE(dat, measurevar="depth_mm", groupvars=c("habitat", "depth"), na.rm=TRUE, .drop=FALSE)
beddat

beddat2 <- beddat %>% mutate(cilow = depth_mm - ci,
                             cihigh = depth_mm + ci)

hist(dat$depth_mm) # very right-skewed
hist(log(dat$depth_mm)) # normal

depthmod <- glmmTMB(log(depth_mm) ~ aspectdepth + 
                      (1|site_unqiue/quadrat), family = "gaussian", data = dat, na.action = "na.omit")

summary(depthmod)

depthmodA <- car::Anova(depthmod) %>% as.data.frame()

write.csv(car::Anova(depthmod), "depthmodA.csv")

depthmod.t <- as_flextable(depthmod) %>%
colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
save_as_docx(depthmod.t, path = "depthmod.t.docx")

pred1 <- emmeans(depthmod, pairwise ~ aspectdepth, type = "response")
# Reef flat thinner bed than all reef slope sites

pred1.1 <- pred1$emmeans %>% as.data.frame()
pred1.1.t2 <- flextable(pred1.1) %>% 
  colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
      save_as_docx(path = "pred1.1.t2.docx")

pred1.2 <- pred1$contrasts %>% as.data.frame()
pred1.2.t <- flextable(pred1.2) %>% 
  colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
  save_as_docx(path = "pred1.2.t.docx")

# Rubble size ----

names(dat)
rublengthdat <- summarySE(dat, measurevar="widest_span_mm", groupvars=c("habitat", "depth"), na.rm=TRUE, .drop=FALSE)
rublengthdat

rublengthdat2 <- rublengthdat %>% mutate(cilow = widest_span_mm - ci,
                             cihigh = widest_span_mm + ci)

hist(dat$widest_span_mm) # very right-skewed
hist(log(dat$widest_span_mm)) # normal

sizemod <- glmmTMB(log(widest_span_mm) ~ aspectdepth +
                     (1|site_unqiue/quadrat), family = "gaussian", data = dat, na.action = "na.omit")

summary(sizemod)

car::Anova(sizemod) # number of pieces is different between aspectdepths
sizemodA <- car::Anova(sizemod) %>% as.data.frame()
sizemodA <- as_flextable(sizemodA)
save_as_docx(sizemodA, path = "sizemodA.docx")

pred1 <- emmeans(sizemod, pairwise ~ aspectdepth, type = "response")
#Size - larger in SE 2-3m compared to reef flat, west 2-3, west 6-7, SE 6-7

pred1.1 <- pred1$emmeans %>% as.data.frame()
pred1.1t <- flextable(pred1.1) %>% 
  colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
  save_as_docx(path = "pred1.1t.docx")

pred1.2 <- pred1$contrasts %>% as.data.frame()
flextable(pred1.2) %>% 
  colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
  set_table_properties(layout = "autofit") %>% 
  save_as_docx(path = "pred1.2t.docx")

# Binding prevalence ----

str(dat)

combmod <- glmmTMB(boundB ~ aspectdepth + morphology_simple + stable + 
                     aspectdepth:morphology_simple +
                     stable:morphology_simple + stable:aspectdepth + 
                     widest_span_mm + widest_span_mm:morphology_simple +
                     widest_span_mm:stable + widest_span_mm:aspectdepth +
                     (1|site_unqiue), family = "binomial", data = dat, na.action = "na.omit") # quadrat cannot be included due to low replication of term combinations

plot(simulateResiduals(fittedModel = combmod)) # good

car::Anova(combmod) 

combmod2 <- update(combmod, ~ . - stable:aspectdepth)
AICc(combmod,combmod2)
car::Anova(combmod2) 
combmod3 <- update(combmod2, ~ . - stable:morphology_simple)
AICc(combmod2,combmod3)
car::Anova(combmod3) 
combmod4 <- update(combmod3, ~ . - aspectdepth:morphology_simple)
AICc(combmod3,combmod4)
car::Anova(combmod4) 
combmod5 <- update(combmod4, ~ . - aspectdepth:widest_span_mm)
AICc(combmod4,combmod5)
car::Anova(combmod5) 
# interaction between size and morphology
# interaction between size and stability
# effect of aspectdepth

plot(simulateResiduals(fittedModel = combmod5)) # good

combmod5A <- car::Anova(combmod5) %>% as.data.frame()
as_flextable(combmod5A) %>% 
  colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
  set_table_properties(layout = "autofit") %>% 
  save_as_docx(path = "combmod5At.docx")

# morphology:size ----

combmod5.grid2 <- with(dat, list(morphology_simple = levels(morphology_simple),
                                 widest_span_mm = c(100, 200, 400)))
combmod5.grid2

range(dat$widest_span_mm,na.rm=TRUE)
mean(dat$widest_span_mm,na.rm=TRUE)

pred1 <- emmeans(combmod5, pairwise ~ morphology_simple | widest_span_mm, 
                 at = combmod5.grid2, type = "response")

pred1 # At small/mean rubble size, there is a difference in morphology,
# but at large sizes, there is no difference
# Branched rubble has a similar prob of binding across sizes, but unbranched and other morphologies
# only have higher probs of binding when rubble pieces are very large

pred1.1 <- pred1$emmeans %>% as.data.frame()
pred1.1t <- flextable(pred1.1) %>% 
  colformat_double(big.mark = ",", digits = 3, na_str = "N/A") %>%
  set_table_properties(layout = "autofit") %>%
  save_as_docx(path = "pred1.1t.docx")

pred1.2 <- pred1$contrasts %>% as.data.frame()
pred1.2t <- flextable(pred1.2) %>% 
  colformat_double(big.mark = ",", digits = 3, na_str = "N/A") %>%
  set_table_properties(layout = "autofit") %>%
  save_as_docx(path = "pred1.2t.docx")

# morphology:size plot

combmod5.grid2 <- with(dat, list(morphology_simple = levels(morphology_simple),
                                 widest_span_mm = seq(min(widest_span_mm, na.rm=TRUE), max(widest_span_mm, na.rm=TRUE), len = 300)))
combmod5.grid2

any(is.na(dat$widest_span_mm))
any(!is.finite(dat$widest_span_mm))

combmod5.preds2 <- emmeans(combmod5, ~ morphology_simple | widest_span_mm, 
                           at=combmod5.grid2, type = 'response') %>% as.data.frame()
combmod5.preds2

pd <- position_dodge(0.5)

colsT3 <- c("Unbranched" = "black", 
            "Branched" = "grey",
            "Other"= "#CDAA7D")

legend_titleT2 <- "Morphology"

dat %>% dplyr::group_by(morphology_simple) %>%
  dplyr::summarise(max = max(widest_span_mm, na.rm=TRUE), 
                   min = min(widest_span_mm, na.rm=TRUE))

chunk1 = combmod5.preds2%>%
  filter((morphology_simple=="Branched" & widest_span_mm < 501)) # highest measured widest span for branched pieces was 501

chunk2 = combmod5.preds2%>%
  filter((morphology_simple=="Other" & widest_span_mm < 451)) # highest measured widest span for stable pieces was 451

chunk3 = combmod5.preds2%>%
  filter((morphology_simple=="Unbranched" & widest_span_mm < 366)) # highest measured widest span for stable pieces was 365

bindGcomb2 <- ggplot(combmod5.preds2, aes(y=prob, x=widest_span_mm, fill = morphology_simple)) +
  geom_ribbon(data=chunk1, aes(x=widest_span_mm, ymax=asymp.UCL, ymin=asymp.LCL),  position = pd, alpha=0.4, show.legend = TRUE) + 
  geom_line(data=chunk1, aes(x=widest_span_mm, y=prob,color=morphology_simple), show.legend = TRUE) +
  geom_ribbon(data=chunk2, aes(x=widest_span_mm, ymax=asymp.UCL, ymin=asymp.LCL),  position = pd, alpha=0.4, show.legend = TRUE) + 
  geom_line(data=chunk2, aes(x=widest_span_mm, y=prob,color=morphology_simple), show.legend = TRUE) +
  geom_ribbon(data=chunk3, aes(x=widest_span_mm, ymax=asymp.UCL, ymin=asymp.LCL),  position = pd, alpha=0.4, show.legend = TRUE) + 
  geom_line(data=chunk3, aes(x=widest_span_mm, y=prob,color=morphology_simple), show.legend = TRUE) +
  labs(x="Widest span (mm)", y = "Probability of rubble being bound") +
  scale_colour_manual(legend_titleT2, values = colsT3, labels=c("Branched", "Other", "Unbranched")) +
  scale_fill_manual(legend_titleT2, values = colsT3, labels=c("Branched", "Other", "Unbranched")) +
  theme_classic() +
  ylim(0,0.9) +
  theme(text = element_text(size=20, face="bold")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position="bottom")

bindGcomb2

ggsave("bindGcomb2.pdf", bindGcomb2)
# Saving 4.42 x 7.96 in image

# stable:size ----

combmod5.grid3 <- with(dat, list(stable = levels(stable),
                                 widest_span_mm = c(100, 150, 180, 190, 200, 400, 500)))
# looking at sizes 10, 15, 18, 19, 20, 40 and 50 cm
combmod5.grid3

pred1 <- emmeans(combmod5, pairwise ~ stable | widest_span_mm, 
                 at = combmod5.grid3, type = "response")

pred1.1 <- pred1$emmeans %>% as.data.frame()
pred1.1t <- flextable(pred1.1) %>% 
  colformat_double(big.mark = ",", digits = 3, na_str = "N/A") %>%
  set_table_properties(layout = "autofit") %>%
  save_as_docx(path = "pred1.1t.docx")

pred1.2 <- pred1$contrasts %>% as.data.frame()
pred1.2t <- flextable(pred1.2) %>% 
  colformat_double(big.mark = ",", digits = 3, na_str = "N/A") %>%
  set_table_properties(layout = "autofit") %>%
  save_as_docx(path = "pred1.2t.docx")

# stable:size plot

combmod5.grid3 <- with(dat, list(stable = levels(stable),
                                 widest_span_mm = seq(min(widest_span_mm, na.rm=TRUE), max(widest_span_mm, na.rm=TRUE), len = 300)))
combmod5.grid3

any(is.na(dat$widest_span_mm))
any(!is.finite(dat$widest_span_mm))

combmod5.preds3 <- emmeans(combmod5, ~ stable | widest_span_mm, 
                           at=combmod5.grid3, type = 'response') %>% as.data.frame
combmod5.preds3

pd <- position_dodge(0.5)

legend_titleT <- "Stability"

colsT2 <- c("Easily mobilised" = "black", 
            "Stable" = "grey")

dat %>% dplyr::group_by(stable) %>%
  dplyr::summarise(max = max(widest_span_mm, na.rm=TRUE), 
                   min = min(widest_span_mm, na.rm=TRUE))
chunk1 = combmod5.preds3%>%
  filter((stable=="Easily mobilised" & widest_span_mm < 386)) # highest measured widest span for easily mobilised pieces was 385

chunk2 = combmod5.preds3%>%
  filter((stable=="Stable" & widest_span_mm < 501)) # highest measured widest span for stable pieces was 500

bindGcomb3 <- ggplot(combmod5.preds3, aes(y=prob, x=widest_span_mm, fill = stable)) +
  geom_ribbon(data=chunk1, aes(x=widest_span_mm, ymax=asymp.UCL, ymin=asymp.LCL),  position = pd, alpha=0.4, show.legend = TRUE) + 
  geom_line(data=chunk1, aes(x=widest_span_mm, y=prob,color=stable), show.legend = TRUE) +
  geom_ribbon(data=chunk2, aes(x=widest_span_mm, ymax=asymp.UCL, ymin=asymp.LCL),  position = pd, alpha=0.4, show.legend = TRUE) + 
  geom_line(data=chunk2, aes(x=widest_span_mm, y=prob,color=stable), show.legend = TRUE) +
  labs(x="Widest span (mm)", y = "Probability of rubble being bound") +
  scale_colour_manual(legend_titleT, values = colsT2, labels=c("Easily mobilised", "Stable")) +
  scale_fill_manual(legend_titleT, values = colsT2, labels=c("Easily mobilised", "Stable")) +
  theme_classic() +
  ylim(0,0.9) +
  theme(text = element_text(size=20, face="bold")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position="bottom")

bindGcomb3

ggsave("bindGcomb3.pdf", bindGcomb3)
# Saving 4.42 x 7.96 in image

# aspectdepth ----

pred1 <- emmeans(combmod5, pairwise ~ aspectdepth, type = "response")

pred1.1 <- pred1$emmeans %>% as.data.frame()
pred1.1t <- flextable(pred1.1) %>% 
  colformat_double(big.mark = ",", digits = 3, na_str = "N/A") %>%
  set_table_properties(layout = "autofit") %>%
  save_as_docx(path = "pred1.1t.docx")

pred1.2 <- pred1$contrasts %>% as.data.frame()
pred1.2t <- flextable(pred1.2) %>% 
  colformat_double(big.mark = ",", digits = 3, na_str = "N/A") %>%
  set_table_properties(layout = "autofit") %>%
  save_as_docx(path = "pred1.2t.docx")

r.squaredGLMM(combmod5) #21%

# aspectdepth plot

combmod5.grid1 <- with(dat, list(aspectdepth = levels(aspectdepth)))
combmod5.grid1

combmod5.preds1 <- emmeans(combmod5, ~ aspectdepth, 
                          at=combmod5.grid1, type = 'response') %>% as.data.frame
combmod5.preds1

pd <- position_dodge(0.5)

colsT <- c("Southeast_2-3m" = "#56B4E9", 
           "Southeast_6-7m" = "#0072B2", 
           "West_2-3m" = "#E69F00", 
           "West_6-7m" = "#D55E00", 
           "Reef flat_2-3m" = "#999999")

#shapesT <- c("Southeast_2-3m" = 19, 
             #"Southeast_6-7m" = 25, 
            # "West_2-3m" = 19, 
            # "West_6-7m" = 25, 
            # "Lagoon_2-3m" = 19)

legend_titleT <- "Habitat"

bindGcomb <- ggplot(combmod5.preds1, aes(y=prob, x=aspectdepth, colour = aspectdepth, fill = aspectdepth)) +
  geom_point(position=pd, size =4) + 
  geom_errorbar(aes(ymin=asymp.LCL, ymax=asymp.UCL), width=.2, position = pd) +
  labs(x="Habitat", y = "Probability of rubble being bound") +
  scale_colour_manual(legend_titleT, values=colsT, labels=c("Reef flat_2-3m", "Southeast_2-3m", "Southeast_6-7m","West_2-3m", "West_6-7m")) +
  scale_fill_manual(legend_titleT, values = colsT, labels=c("Reef flat_2-3m", "Southeast_2-3m", "Southeast_6-7m","West_2-3m", "West_6-7m")) +
  #scale_shape_manual(legend_titleT, values = shapesT, labels=c("Lag_Shal", "Shelt_Shal", "Shelt_Deep","Exp_Shal", "Exp_Deep")) +
  theme_classic() +
  ylim(0,0.6) +
  theme(text = element_text(size=20, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position="bottom")

bindGcomb

ggsave("bindGcomb.pdf", bindGcomb)
# Saving 4.42 x 7.96 in image

# Stability between habitats ----

names(dat)

stabmod <- glmmTMB(stable ~ aspectdepth + morphology_simple + widest_span_mm +
                     aspectdepth:morphology_simple +
                     widest_span_mm:morphology_simple + 
                     widest_span_mm:aspectdepth +
                     (1|site_unqiue/quadrat), family = "binomial", data = dat, na.action = "na.omit")

car::Anova(stabmod)

stabmod2 <- update(stabmod, ~ . - aspectdepth:morphology_simple)

car::Anova(stabmod2)
AICc(stabmod, stabmod2)

stabmod3 <- update(stabmod2, ~ . - morphology_simple:widest_span_mm)

car::Anova(stabmod3)

write.csv(car::Anova(stabmod3), "car::Anova(stabmod2).csv")

# Morphology

preds1 <- emmeans(stabmod3, pairwise ~ morphology_simple, type = "response")

write.csv(preds1$emmeans, "morphologyemmeans.csv")

write.csv(preds1$contrasts, "morphologycontrasts.csv")

# Aspectdepth and widest span

stabmod2.grid <- with(dat, list(widest_span_mm = c(100, 200),
                                aspectdepth = levels(aspectdepth)))

preds2 <- emmeans(stabmod3, pairwise ~ widest_span_mm | aspectdepth, 
        at = stabmod2.grid, type = "response") # stability increases with size, as expected.

write.csv(preds2$emmeans, "sizeemmeans.csv")

write.csv(preds2$contrasts, "sizeemmeans2.csv")

preds3 <- emmeans(stabmod3, pairwise ~  aspectdepth | widest_span_mm, 
                  at = stabmod2.grid, type = "response") # stability increases with size, as expected.

write.csv(preds3$emmeans, "sizeemmeans3.csv")

write.csv(preds3$contrasts, "sizeemmeans3_2.csv")

# Community composition ----

datY <- filter(dat, boundB == "1") %>% droplevels() # only looking at pairs that are bound

cols1 <- c("habitat", "site", "depth", "trip", "remove_for_size", "quadrat","square",
           "morphology_simple", "stable", "bound", "binder_1", "binder_2", "binder_3")

datY <- datY %>%  mutate_at(cols1, factor)

levels(datY$binder_1)
levels(datY$binder_2)
levels(datY$binder_3)

datY$binder_1 <- factor(datY$binder_1,
                         levels = c("anemone" , "ascidian", "brozoan" ,  "cca",   "coral", "macroalgae",
                                    "old_cemented" ,"serpulid" ,   
                                     "sponge"   ,    "turf"      ,   "unknown"), 
                         labels = c("anemone" , "ascidian", "bryozoan" ,  "cca",   "coral", "macroalgae",
                                    "old_cemented" ,"serpulid" ,   
                                    "sponge"   ,    "turf"      ,   "unknown"))

# Data exploration

ggplot(datY, aes(aspectdepth)) +
  geom_bar(aes(fill=binder_1), position = "fill")

ggplot(datY, aes(aspectdepth)) +
  geom_bar(aes(fill=binder_1), position = "stack")

# Put in long form

datY_long <- gather(datY, key = "bindernumber", value = "type", 
                    binder_1, binder_2, binder_3, na.rm=TRUE)

View(datY_long)

datY_long$type <- as.factor(datY_long$type)
levels(datY_long$type)

# How many of each?

stats <- datY_long %>% group_by(type) %>% 
  dplyr::summarise(N=n(), Percent = (N/497)*100) %>% as.data.frame()
stats

datY_long %>% group_by(bindernumber) %>% 
  dplyr::summarise(N=n()) %>% as.data.frame()

sum(stats$N)
# 497 dominant binders in total. 355 pieces of bound rubble
# 237 pieces were bound by one binder type only, 94 by 2 types, and 24 by 3 dominant binders.

# Across different habitats and depths, sponge and coralline algae were the most prolific binders.

# 46.2% of bound pairs had sponge as a dominant binder (at least one sponge bind, individual binds were not counted)
# 31.5% of bound pairs had CCA as a dominant binder
# 16.6% of bound pairs had at least one bind that appeared to be quite old and could have been 
# cemented following CCA binding, or potentially, the coral pieces had anastomised when living (Cemented category)
# 17.5% of bound pairs had at least one macroalgae bind (common macroalgae binders were Dictyota, and Tydemania)
# 6.8% of bound pairs had at least one coral bind 

stats2 <- dat %>% group_by(aspectdepth, bound) %>% 
  dplyr::summarise(N=n()) %>% as.data.frame()
stats2
#habitat   depth total
#1 Lagoon    2-3m     28 bound rubble pieces. 42 dominant binders.
#2 Southeast 2-3m     64 bound rubble piece. 80 dominant binders.
#3 Southeast 6-7m     65 bound rubble pieces. 85 dominant binders.
#4 West      2-3m    82 bound rubble pieces. 127 dominant binders.
#5 West      6-7m    116 bound rubble pieces. 163 dominant binders

stats4 <- datY_long %>% group_by(aspectdepth, bindernumber, type) %>% 
  dplyr::summarise(N=n()) %>% as.data.frame()
stats4

stats3 <- datY_long %>% group_by(aspectdepth, type) %>% 
  dplyr::summarise(N=n()) %>% as.data.frame()
stats3

12/42 #29% of dominant binds in the reef flat were attributed to CCA. 12/28 bound rubble pieces = 43%
9/42 #21% of binds were solid. 9 rubble pieces out of 28 had solid binds = 32%
# In the reef flat, the most represented binder type was coralline algae
8/42 #19% sponge 8/28 rubble pieces = 29%

# Reef slope.
# SE 2-3
27/80 #34% sponge. 27/64 = 42%
14/80 #18% of binders cca ; 14 rubble pieces. 14/64 = 22%
15/80 #19 macroalgae ; 15/64 = 23%
11/80 #14 solid 11/64 = 17%
9/80 #11 ascidians ; 9/64 = 14%

# SE 6-7
38/85 #45% sponge bind. 38/65 = 58%
15/85 #18 cca ; 15/65 = 23%
10/85 #12 macroalgae ; 10/65 = 15%
8/85 #9 solid 8/65 = 12%
4/85 #5 ascidians 4/65 = 6%

# W 2-3
34/127 # 27% sponge. ; 34/82 41%
20/127 # 16% coral. ; 20/82 24%
#*Looking at the raw data, 19 rubble pairs had at least 1 coral dominant binder. 19/82 bound pieces = 23%**
23/127 #18 cca ; 23/82 28%
16/127 #13 macroalgae ; 16/82 20%
14/127 #11 solid ; 14/82 17%
9/127 #7 ascidians ; 9/82 11%

# W 6-7
57/163 # 35% sponge.; 57/116 49%
48/163 # 30 cca ; 48/116 41%
21/163 # 13 macroalgae ; 21/116 18%
17/163 #10 solid ; 17/116 15%
8/163 #5 ascidians 8/116 7%

# On the reef slope, the most represented binders were sponges.
# Bound pieces by coral were found in SE 6-7m sites and both depths on the western side of the island,
# being most prolific at W 2-3m:
# Here, 16% of binds were from corals (predominantly Porites rus and Leptoseris)

View(dat)

# Binders plot -----

levels(datY_long$type)

datY_long$type <- factor(datY_long$type,
                         levels = c('turf', 'macroalgae', 'cca',
                                    'sponge', 'ascidian', 'bryozoan',
                                    'anemone', 'coral', 'serpulid',
                                    'unknown', 'old_cemented'), 
                         labels = c('Turf algae', 'Macroalgae', 'Coralline algae',
                                    'Sponge', 'Ascidian', 'Bryozoan',
                                    'Anemone', 'Hard coral', 'Serpulidae',
                                    'Other', 'Solid'))

typecols <- c("Turf algae" = "#B4EEB4",
              "Macroalgae" = "#79B07A",
              "Coralline algae" =  "#EEB4B4",
              "Sponge" = "#F0E68C",     
              "Ascidian" = "#C799C6",
              "Bryozoan" =   "lightblue2",
              "Anemone" = "#FA8072",
              "Hard coral" = "peachpuff3",
              "Serpulidae" = "#EEE8CD",
              "Other" =  "#C7C7C7",
              "Solid" = "#919DBD")

bindcols <- c("Turf algae" = "darkseagreen1",
              "Macroalgae" = "aquamarine3",
              "Coralline algae" =  "lightpink",        
              "Sponge" = "khaki1",     
              "Ascidian" = "mediumpurple1",
              "Bryozoan" =   "lightblue2",
              "Anemone" = "indianred1",
              "Hard coral" = "peachpuff3",
              "Serpulidae" = "ivory",
              "Other" =  "grey",
              "Solid" = "#919DBD")

typebar <- ggplot(datY_long, aes(x=aspectdepth)) +
  geom_bar(aes(fill = type), position="fill", colour="black") + 
  labs(x="Habitat", y = "Proportion binder type") + 
  scale_fill_manual(name = "Binder type", values= typecols) +
  scale_x_discrete(labels = c("Reef flat 2-3m", "Sheltered 2-3m", "Sheltered 6-7m",
                              "Exposed 2-3m", "Exposed 6-7m")) +
  # scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
  #                   breaks = c(0.0,0.2,0.4,0.6,0.8,1.0),
  #                  limit = c(-0.05,1.0)) +
  theme_classic() +
  #ggtitle("Binder composition of bound rubble in in-situ rubble beds") +
  theme(text = element_text(face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "bottom") +
  theme(axis.text.x.bottom = element_text(size=20,vjust = 0.6,hjust = 0.5,angle = 45))

typebar

ggsave("typebar3.pdf", typebar)
# Saving 11.4 x 7.96 in image

# Plot with total count

ggplot(datY_long, aes(x=aspectdepth)) +
  geom_bar(aes(fill = type), position="stack") + 
  labs(x="Habitat", y = "Proportion binder type") + 
  scale_fill_manual(name = "Binder type", values= typecols) +
  scale_x_discrete(labels = c("Reef flat 2-3m", "Sheltered 2-3m", "Sheltered 6-7m",
                              "Exposed 2-3m", "Exposed 6-7m"))


# PERMANOVA ----

# Formatting the dataframe prior to being able to run a permanova on a distance matrix

# Put datlong into wide form
# First need to group by and get a count per organism
# Will not include the Not bounds because they are causing problems in the distance matrix analysis
# if they are left in

datY_long2 <- gather(datY, key = "bindernumber", value = "type", 
                    binder_1, binder_2, binder_3)

# Format types
datY_long2$type <- as.factor(datY_long2$type)
levels(datY_long2$type)

datY_long2$type <- factor(datY_long2$type,
                         levels = c('turf', 'macroalgae', 'cca',
                                    'sponge', 'ascidian', 'bryozoan',
                                    'anemone', 'coral', 'serpulid',
                                    'unknown', 'old_cemented'), 
                         labels = c('Turf_algae', 'Macroalgae', 'Coralline_algae',
                                    'Sponge', 'Ascidian', 'Bryozoan',
                                    'Anemone', 'Hard_coral', 'Serpulid_worm',
                                    'Other', 'Solid'))
levels(datY_long2$type)
view(datY_long2)

# Now we want the total per type, before we can put it back in wide form 
# We want wide form with categories of binders as columns instead of binder_1, binder_2, binder_3 columns
datY_long3 <- datY_long2 %>% dplyr::group_by(habitat, site, trip, quadrat2,
                                     piece, aspectdepth, type) %>%
                    dplyr::summarise(typecount = n())

view(datY_long3)
levels(datY_long3$type)

# Remove cases where type = NA, these are the ones where there might have been a binder_1, but
# no binder_2 or binder_3, etc.

datY_long4 <- datY_long3 %>%
  filter(!is.na(type))

# Now put into wide form again
datY_wide2 <- datY_long4 %>% spread(key = type, value = typecount)
view(datY_wide2)

# Change all the NAs to 0s
datY_wide2[is.na(datY_wide2)] <- 0

view(datY_wide2)

# Make community and meta matrices

datY_wide2com <- datY_wide2 %>% ungroup() %>% #Makes the community matrix (the response variables)
  dplyr::select(Turf_algae, Macroalgae, Coralline_algae, Sponge,         
                Ascidian, Bryozoan, Anemone, Hard_coral, Serpulid_worm,  
                Other, Solid)

# Will change matrix to just presence or absence because there are only a few cases where numbers are above 1

datY_wide2com <- datY_wide2com %>% mutate_all( ~ifelse(.>0, 1, 0))
View(datY_wide2com)

datY_wide2met <- datY_wide2 %>% #Makes a meta matrix (the explanatory variables)
  dplyr::select(habitat, site, trip, quadrat2, piece, aspectdepth)

names(datY_wide2met)

datY_wide2com.dist <-vegdist(datY_wide2com, method = "jaccard")

# Homogeneity of variance check
b=betadisper(datY_wide2com.dist,group=datY_wide2met[,"aspectdepth"]) #checking homogeneity of variance of your factors. Distance matrix is first argument following by factor of your meta matrix (the factors)
anova(b)

TukeyHSD(b)

# Permanova is robust to violations of homogeneity of variance IF the design is balanced, which this is.

# Fit the model

comm.adonis <-adonis2(datY_wide2com.dist ~ aspectdepth, 
                         data = datY_wide2, method = "jaccard", permutations = 999) %>% as.data.frame()

comm.adonis

write.csv(comm.adonis, "comm.adonis.csv")

perm1 <- flextable(comm.adonis %>% as.data.frame()) %>% 
  colformat_double(big.mark = ",", digits = 2, na_str = "N/A") %>%
  save_as_docx(path = "perm1.docx")

compsN <- pairwiseAdonis::pairwise.adonis(datY_wide2com[,1:11],datY_wide2met$aspectdepth) 

compsN

write.csv(compsN, "compsN.csv")

# Simper

(simN <- with(datY_wide2met, simper(datY_wide2com, aspectdepth, permutations = 999)))
summary(simN)


# nMDS ----

# nMDS of the natural binder community (presence / absence of the top 3 binders)

bindcomnmds <- metaMDS(datY_wide2com, distance = "jaccard", k=2,trymax=200, autotransform=TRUE,expand=FALSE, plot=FALSE)

bindcomnmds$stress # The stress is 0.04. A value below 0.2 is acceptable, stress is 1-R2 value.

stressplot(bindcomnmds) # Shows that our non-metric fit R2 is 0.99. There is just not a lot of data!!

bindcomnmds.sites.scores<- as.data.frame(scores(bindcomnmds, display = 'sites'))
bindcomnmds.sites.scores<- data.frame(bindcomnmds.sites.scores, datY_wide2met)
bindcomnmds.species.scores<- as.data.frame(scores(bindcomnmds, display = 'species'))
bindcomnmds.species.scores$Binders<- row.names(bindcomnmds.species.scores)

head(bindcomnmds.sites.scores)
head(bindcomnmds.species.scores)

bindcomnmds.species.scores$Binders <- as.factor(bindcomnmds.species.scores$Binders)
levels(bindcomnmds.species.scores$Binders)

bindcomnmds.species.scores$Binders <- factor(bindcomnmds.species.scores$Binders,
                                             levels = c("Turf_algae",
                                                        "Macroalgae",
                                                        "Coralline_algae",                
                                                        "Sponge",             
                                                        "Ascidian",
                                                        "Bryozoan",
                                                        "Anemone",
                                                        "Hard_coral", 
                                                        "Serpulid_worm", 
                                                        "Other",
                                                        "Solid"),    
                                             labels = c("Turf_algae",
                                                        "Macroalgae",
                                                        "Coralline_algae",                
                                                        "Sponge",             
                                                        "Ascidian",
                                                        "Bryozoan",
                                                        "Anemone",
                                                        "Hard_coral", 
                                                        "Serpulid_worm", 
                                                        "Other",
                                                        "Solid"))

bindcomnmds.sites.scores$aspectdepth <- as.factor(bindcomnmds.sites.scores$aspectdepth)

levels(bindcomnmds.sites.scores$aspectdepth)

bindcomnmds.sites.scores$aspectdepth <- 
  factor(bindcomnmds.sites.scores$aspectdepth, 
         levels = c("Reef flat_2-3m", "Southeast_2-3m", "Southeast_6-7m", "West_2-3m", "West_6-7m" ),
         labels = c("RF", "SheltShal", "SheltDeep", "ExpShal", "ExpDeep"))

colsT <- c("SheltShal" = "#56B4E9",
           "SheltDeep" = "#0072B2",
         "ExpShal" = "#E69F00",
           "ExpDeep" = "#D55E00",
           "RF" = "#999999")

shapesT <- c("SheltShal" = 19,
             "SheltDeep" = 25,
             "ExpShal" = 19,
             "ExpDeep" = 25,
             "RF" = 19)

centroidsA <- aggregate(cbind(NMDS1, NMDS2) ~ aspectdepth, data = bindcomnmds.sites.scores, FUN = mean)
centroidsA <- data.frame(centroidsA) # pulling out just the averages which is what I will plot
legend_title <- "Habitat & Month"

library(ggrepel)

# Habitat points ----

plotbindcomm <-ggplot()+
  geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
  geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')+
  #geom_segment(data=topsideshallnmds.species.scores, aes(y=0, x=0, yend=NMDS2, xend=NMDS1), colour = "grey", arrow=arrow(length=unit(0.3, 'lines')))+ #adding the binding organism drivers
  #geom_text(data=topsideshallnmds.species.scores, aes(y=NMDS2, x=NMDS1, label = Binders, hjust="inward",vjust="inward"), show.legend=FALSE, colour = "grey") +
  geom_point(data=centroidsA, aes(y=NMDS2, x=NMDS1, color=aspectdepth, shape = aspectdepth, fill = aspectdepth), size = 3)+
  scale_colour_manual(legend_title, values=colsT) +
  scale_fill_manual(legend_title, values = colsT) +
  scale_shape_manual(legend_title, values = shapesT) +
  geom_text_repel(data=centroidsA, aes(y=NMDS2, x=NMDS1, label = aspectdepth,
                                       hjust="inward",vjust="inward", color=aspectdepth), show.legend=FALSE) +
  #geom_segment(data=bindcomnmds.species.scores, aes(y=0, x=0, yend=NMDS2, xend=NMDS1), 
              # colour = "grey", arrow=arrow(length=unit(0.3, 'lines'))) + 
  #geom_text_repel(data=bindcomnmds.species.scores, aes(y=NMDS2, x=NMDS1, 
                                                      # hjust="inward",vjust="inward",
                                                       #label = Binders), colour = "gray27") + 
  theme_classic() +
  ylab("NMDS2") +
  xlab("NMDS1") +
  ggtitle("In situ Binder Community") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_vline(xintercept=0, linetype="dashed") + 
  geom_hline(yintercept=0, linetype="dashed") +
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "none") 
plotbindcomm

# Binder segments plot ----

plotbindseg <-ggplot()+
  geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
  geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')+
  geom_segment(data=bindcomnmds.species.scores, aes(y=0, x=0, yend=NMDS2, xend=NMDS1), 
               colour = "grey", arrow=arrow(length=unit(0.3, 'lines'))) + 
  geom_text_repel(data=bindcomnmds.species.scores, aes(y=NMDS2, x=NMDS1, 
                                                       hjust="inward",vjust="inward",
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
plotbindseg

grid.arrange (plotbindcomm, plotbindseg, nrow = 2)


# Graphical abstract ----

# Prop of rubble bound

dat.sum2 <- dat %>%
  group_by(habitat, aspectdepth) %>%
  dplyr::summarise(meanB = mean(boundB, na.rm=TRUE),
                   seB = sd(boundB, na.rm=TRUE)/sqrt(n()),
                   N=n())

# Reef flat is 7% bound
# Southeast is 20% bound, so 1/5 rubbles black
# West is 30% so 2/5 rubbles black

# Number of binds per exposure and month

# Not possible to estimate for the surveys, as only the dominant binders were recorded, not the 
# number of binds, so, we have number of binder groups (up to 3 though)

# Prop of binds of each binder type (on average)

names(datY_long)
stats4 <- datY_long %>% group_by(habitat, type) %>% 
  dplyr::summarise(N=n()) %>% as.data.frame()
stats4

stats5 <- datY_long %>% dplyr::group_by(habitat) %>% 
  dplyr::summarise(N=n())

# None of the groups make up 50% or more on their own


