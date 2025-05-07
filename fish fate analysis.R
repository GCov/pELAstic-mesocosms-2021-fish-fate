#### Load Packages ####

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(randomForest)
library(caret)

theme1 <-
  theme_bw() +
  theme(
    panel.spacing = unit(0.5, "cm"),
    text = element_text(size = 9,
                        family = "serif"),
    axis.text = element_text(size = 9),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing = unit(0, "cm")
  )


#### Load Data and Prepare ####

# First load in the fish data

fish <- read.csv("fish_data.csv",
                 header = TRUE,
                 stringsAsFactors = TRUE)

str(fish)

# Add treatments and change name of organ from fillet to muscle

treatments <- 
  data.frame(corral = as.factor(c("B", "C", "D", "E", "F", "G", "H", "I")),
             nominal_MPs = as.numeric(c(0, 414, 29240, 100, 6, 7071, 0, 1710)))

fish2 <- 
  fish %>% 
  left_join(treatments,
            by = "corral")

levels(fish2$organ) <- c("Muscle", "Gill", "GIT", "Liver")

# Add biometrics data

biometrics <- read.csv("fish_biometrics.csv",
                       header = TRUE,
                       stringsAsFactors = TRUE)

str(biometrics)


fish2 <-
  fish2 %>% 
  left_join(biometrics,
             by = "fish_ID")

# Those data are ready to go
# Load in blanks and spike and recovery data
# We'll deal with these data later


blanks <- read.csv("blanks.csv", 
                   header = TRUE,
                   stringsAsFactors = TRUE)

str(blanks)

SaR <- read.csv("spike_and_recovery.csv",
                header = TRUE,
                stringsAsFactors = TRUE)

str(SaR)

#### Spectroscopy Correction ####

# Try correcting the data using a random forest model to predict whether
# particles should be 'real' ELA particles or not

##### Set up the data #####

# First separate out just the particles and look at them in a few ways

fish_particles <-
  fish2 %>% 
  filter(!is.na(particle_number))

ggplot(fish_particles, 
       aes(x = spectroscopy_method, 
           fill = organ)) + 
  geom_bar() + 
  facet_wrap(~ count_ID, scales = "free_y") + 
  labs(x = "Spectroscopy Method", y = "Number of particles")

# Now separate out just particles with spectroscopy IDs

spectroscopy_particles <-
  fish_particles %>% 
  filter(!is.na(spectroscopy_method))

# Explore the data a little bit

summary(spectroscopy_particles$count_ID)
summary(spectroscopy_particles$spectroscopy_method)  # decent data for each method
summary(spectroscopy_particles$spectroscopy_ID)

spectroscopy_summary<- 
  spectroscopy_particles %>% 
  group_by(count_ID, spectroscopy_ID, organ) %>% 
  summarize(count = length(spectroscopy_ID))

ggplot(spectroscopy_particles, 
       aes(x = count_ID, 
           fill = spectroscopy_ID)) + 
  geom_bar() + 
  facet_grid(spectroscopy_method ~ organ, scales = "free_y") + 
  labs(x = "Polymer Identified by Colour", y = "Number of particles") +
  scale_fill_viridis_d(option = "turbo",
                       name = "Spectroscopy ID")

# Make a new column for category

spectroscopy_particles <-
  spectroscopy_particles %>% 
  mutate(cat = as.factor(ifelse(spectroscopy_ID == "PE", 
                                "PE", 
                                ifelse(spectroscopy_ID == "PET",
                                       "PET",
                                       ifelse(spectroscopy_ID == "PS", 
                                              "PS", 
                                              "Contamination")))))

# Plot

ggplot(spectroscopy_particles, 
       aes(x = count_ID,
           fill = cat)) +
  geom_bar() +
  labs(x = "Microscopy ID",
       y = "Count") +
  scale_fill_manual(name = "Spectroscopy Class",
                    values = c("grey", "yellow", "blue", "pink")) +
  facet_wrap(~ organ, scales = "free_y") +
  theme_bw()

# Split into 70% training 70% test

set.seed(42345)

train <-
  spectroscopy_particles %>% 
  slice_sample(prop = 0.7, by = c(organ, cat))

test <-
  spectroscopy_particles %>% 
  anti_join(train,
            by = c("fish_ID", "particle_number"))

# Make sure that this worked

ggplot(train, aes(x = cat,
                  y = particle_number,
                  fill = organ)) +
  geom_col() +
  theme_bw()

ggplot(test, aes(x = cat,
                  y = particle_number,
                  fill = organ)) +
  geom_col() +
  theme_bw()

# Seems to have worked pretty well (representative sampling)

##### Random Forest Model ####

set.seed(42354)

rf <- randomForest(cat ~ organ + count_ID + fish_ID,
                   data = train,
                   proximity = TRUE,
                   ntree = 1501)

rf  # OOB error rate = 9.64%%

plot(rf)

rf2 <- randomForest(cat ~ organ + count_ID + fish_ID,
                   data = train,
                   proximity = TRUE,
                   ntree = 251)

rf2  # OOB error rate = 10.04%%

train$predict <- predict(rf2)

test$predict <- predict(rf2, test)

CMtrain <-
  as.data.frame(confusionMatrix(train$predict,
                                train$cat)$table)

ggplot(CMtrain, 
       aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = 'red') +
  scale_fill_viridis_c(option = "mako",
                       direction = -1) +
  labs(x = "Data",
       y = "Prediction") +
  theme_minimal()

# This looks pretty good!!!

CMtest <-
  as.data.frame(confusionMatrix(test$predict,
                                test$cat)$table)

ggplot(CMtest, 
       aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = 'red') +
  scale_fill_viridis_c(option = "mako",
                       direction = -1) +
  labs(x = "Data",
       y = "Prediction") +
  theme_minimal()

# This also looks quite good!

varImpPlot(rf2)  # count_ID most important

`# This seems OK so refit the model to the full data set

set.seed(424)

rf3 <- randomForest(cat ~ count_ID + organ,
                    data = spectroscopy_particles,
                    proximity = TRUE,
                    ntree = 1501)

rf3  # OOB error rate = 9.32%

plot(rf3)

rf4 <- randomForest(cat ~ count_ID + organ,
                    data = train,
                    proximity = TRUE,
                    ntree = 1001)

rf4  # OOB error rate = 9.24%

###### Adjust Particles with no spectroscopy ####

no_spectroscopy_particles <-
  fish_particles %>% 
  filter(is.na(spectroscopy_method))

no_spectroscopy_particles$predicted_class <-
  predict(rf4, no_spectroscopy_particles)

summary(no_spectroscopy_particles$predicted_class)
# assigned 67 particles as contamination

# Plot how the assignments match up with microscopy IDs

ggplot(no_spectroscopy_particles, 
       aes(x = count_ID,
           fill = predicted_class)) +
  geom_bar() +
  labs(x = "Microscopy ID",
       y = "Count") +
  scale_fill_manual(name = "Predicted Class",
                    values = c("grey", "yellow", "blue", "pink")) +
  facet_wrap(~ organ, scales = "free_y") +
  theme_bw()

# So basically the only 'correction' this has made is to remove all PS from the
# fillet samples that didn't have spectroscopy done.

fish3 <-
  fish2 %>% 
  left_join(no_spectroscopy_particles) %>%
  left_join(spectroscopy_particles) %>% 
  mutate(polymer = coalesce(cat, predicted_class))

# Plot all particles

ggplot(fish3 %>% filter(!is.na(count_ID)), 
       aes(x = count_ID,
           fill = polymer)) +
  geom_bar() +
  labs(x = "Microscopy ID",
       y = "Count") +
  scale_fill_manual(name = "Assigned Polymer",
                    values = c("grey", "yellow", "blue", "pink")) +
  facet_wrap(~ organ, scales = "free_y") +
  theme_bw()

# Now I can work with the data as if we know the polymer class for each
# individual particles

##### Plot the Adjusted Data ####

ggplot(fish3 %>% filter(!is.na(polymer)),
       aes(x = organ, 
           fill = polymer)) +
  geom_bar() +
  facet_wrap(~ nominal_MPs, scales = "free_y") +
  scale_fill_manual(name = "Assigned Polymer",
                    values = c("grey", "yellow", "blue", "pink")) +
  theme_bw()

ggplot(fish3 %>% filter(!is.na(polymer)),
       aes(x = organ,
           fill = polymer)) +
  geom_bar() +
  facet_wrap(~ fish_ID, scales = "free_y") +
  scale_fill_manual(name = "Assigned Polymer",
                    values = c("grey", "yellow", "blue", "pink")) +
  theme_bw()

#### Summarize New Data ####

# First remove gill data and contamination particles
# Then summarize by particle count by polymer

fish4 <-
  fish3 %>% 
  filter(organ != "Gill") %>% 
  mutate(particle = ifelse(polymer == "PET" | 
                             polymer == "PE" |
                             polymer == "PS",
                           1,
                           0))

fish4_particles <-
  fish4 %>% 
  filter(!is.na(particle))

# Need to make a new data set with every polymer for every sample then
# populate it

fish_new <-
  expand.grid(fish_ID = unique(fish4$fish_ID),
              organ = unique(fish4$organ),
              polymer = unique(fish4$polymer)[-c(1,4)])

# Add in sample data

fish_data <-
  fish4 %>% 
  group_by(corral,
           fish_ID,
           organ,
           nominal_MPs,
           total_length,
           weight,
           GIT_weight,
           fillet_weight,
           liver_weight,
           gill_weight) %>% 
  summarize(dummy = length(sample_ID),
            .groups = "drop")

# Combine into new data

fish_new <-
  fish_new %>% 
  left_join(fish_data, 
            by = c("fish_ID", "organ"))

# Summarise counts of fish particles

fish_particle_summary <-
  fish4_particles %>% 
  group_by(fish_ID,
           organ,
           polymer) %>% 
  summarize(count = sum(particle),
            .groups = "drop") %>% 
  ungroup()

fish_full_summary <-
  fish_new %>% 
  left_join(fish_particle_summary,
            by = c("fish_ID",
                   "organ",
                   "polymer")) %>% 
  mutate(count = replace_na(count, 0))

# Plot

ggplot(fish_full_summary,
       aes(x = as.factor(as.character(nominal_MPs)),
           y = count,
           fill = polymer)) +
  geom_boxplot() +
  geom_point(shape = 21) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 100, 1000, 10000)) +
  facet_grid(organ ~ polymer,
             scales = "free_y") +
  labs(y = "MPs per ind",
       x = "Nominal Exposure Concentration MPs") +
  scale_fill_manual(values = c("yellow", "blue", "pink"),
                    name = "Polymer") +
  theme_bw()

# Summary stats

fish_totals <-
  fish_full_summary %>% 
  group_by(fish_ID, organ, nominal_MPs) %>% 
  summarize(total_count = sum(count),
            .groups = "drop")

fish_means <-
  fish_totals %>% 
  group_by(organ, nominal_MPs) %>% 
  summarize(mean = mean(total_count),
            sd = sd(total_count),
            .groups = "drop")

fish_means

# Plot

ggplot(fish_totals,
       aes(x = nominal_MPs,
           y = total_count,
           group = nominal_MPs, 
           fill = as.factor(nominal_MPs))) +
  geom_boxplot() +
  geom_jitter(shape = 21,
              height = 0) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  facet_wrap(~ organ,
             scales = "free_y") +
  labs(y = "MPs per ind",
       x = "Nominal Exposure Concentration MPs") +
  scale_fill_viridis_d(option = "inferno",
                       direction = -1,
                       name = "") +
  theme_bw()

#### Blank and Recovery Correction ####

##### Explore the blanks data #####

blanks

# 12 blanks; 1 PE particles in 1 blank sample

# Summarize the blanks data

blanks2 <-
  blanks %>% 
  mutate(particle = ifelse(is.na(particle_number), 
                           0,
                           1),
         polymer = ID)

blanks_new <-
  expand.grid(sample_ID = unique(blanks$sample_ID),
              polymer = c("PE", "PS", "PET"))

blanks_summary <- 
  blanks_new %>% 
  left_join(blanks2,
             by = c("sample_ID",
                    "polymer")) %>% 
  mutate(count = replace_na(particle, 0))

# Calculate means

blank_means <-
  blanks_summary %>% 
  group_by(polymer) %>% 
  summarize(mean_blanks = mean(count),
            sd_blanks = sd(count))

blank_means

##### Explore the spike and recovery data #####

# Note that I modified the raw data a bit to get it into a format I liked
# This included deleting on PS particle that was beyond the 10 particles put in
# Note this in the discussion though - recovery was 11/10 for SR3 PS

head(SaR)

# Summarize by sample/polymer

SaR_summary <-
  SaR %>% 
  group_by(sample_ID,
           polymer) %>% 
  summarize(recovery = length(recovery[recovery == "yes"]),
            .groups = "drop")

# Plot this

ggplot(SaR_summary, 
       aes(x = sample_ID,
           y = recovery,
           fill = polymer)) +
  geom_col() +
  facet_wrap(~ polymer) +
  scale_fill_manual(values = c("yellow", "blue", "pink")) +
  theme_bw()

# Calculate means

SaR_means <-
  SaR_summary %>% 
  group_by(polymer) %>% 
  summarize(mean_prob_recovery = mean(recovery / 10),
            sd_prob_recovery = sd(recovery / 10))

#### Modeling ####

# See how organ concentration relate to nominal exposure concentrations

# I need to generate a new column in the data for weight of the specific organ
# for each row

fish_full_summary$organ_weight <- 0

for(i in 1:nrow(fish_full_summary)) {
  fish_full_summary$organ_weight[i] <-
    ifelse(fish_full_summary$organ[i] == "Fillet",
           fish_full_summary$fillet_weight[i],
           ifelse(fish_full_summary$organ[i] == "GIT",
                  fish_full_summary$GIT_weight[i],
                  fish_full_summary$liver_weight[i]))
}

fish_full_summary <-
  fish_full_summary %>% 
  mutate(log_organ_weight = log(organ_weight))

# grid for prediction
reference_grid <-
  expand.grid(nominal_MPs = unique(fish_full_summary$nominal_MPs),
              polymer = unique(fish_full_summary$polymer))

##### GLMM with addition/subtraction #####

## Modified to use separate models for each organ

###### Corrected data #####

fish_full_summary2 <-
  fish_full_summary %>% 
  left_join(SaR_means,
            by = "polymer") %>% 
  left_join(blank_means,
            by = "polymer") %>% 
  mutate(adjusted_count = floor((count / mean_prob_recovery) - mean_blanks))

fish_full_summary2$adjusted_count[fish_full_summary2$adjusted_count < 0] <- 0

# Compare original and adjusted data

ggplot(fish_full_summary2) +
  geom_point(aes(x = count, y = adjusted_count)) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  facet_wrap(polymer ~ organ, scales = "free")

# Separate data by organ

GITs <-
  fish_full_summary2 %>% 
  filter(organ == "GIT")

livers <- 
  fish_full_summary2 %>% 
  filter(organ == "Liver")

muscle <-
  fish_full_summary2 %>% 
  filter(organ == "Muscle")

# Now fit the models

###### GITs ####

# Poisson, no interactions

GIT_mod1 <-
  glmmTMB(adjusted_count ~ 0 + polymer + log(nominal_MPs + 6) + 
            log(total_length) + (1 | corral),
                family = poisson(link = "log"),
                data = GITs)

plot(simulateResiduals(GIT_mod1))  # fail

# Poisson, 2-way interaction

GIT_mod2 <-
  glmmTMB(adjusted_count ~ 0 + polymer * log(nominal_MPs + 6) + 
            log(total_length) + (1 | corral),
          family = poisson(link = "log"),
          data = GITs)

plot(simulateResiduals(GIT_mod2))  # fail

# Negative binomial, no interactions

GIT_mod3 <-
  glmmTMB(adjusted_count ~ 0 + polymer + log(nominal_MPs + 6) + 
            log(total_length) + (1 | corral),
          family = nbinom2(link = "log"),
          data = GITs)

plot(simulateResiduals(GIT_mod3))  # fail

# Negative binomial, 2-way interaction

GIT_mod4 <-
  glmmTMB(adjusted_count ~ 0 + polymer * log(nominal_MPs + 6) + 
            log(total_length) + (1 | corral),
          family = nbinom2(link = "log"),
          data = GITs)

plot(simulateResiduals(GIT_mod4))  # fail

# Negative binomial, 3-way interaction

GIT_mod5 <-
  glmmTMB(adjusted_count ~ 0 + polymer * log(nominal_MPs + 6) * 
            log(total_length) + (1 | corral),
          family = nbinom2(link = "log"),
          data = GITs)

plot(simulateResiduals(GIT_mod5))  # fail

# ZIP, 2-way interaction

GIT_mod6 <-
  glmmTMB(adjusted_count ~ 0 + polymer * log(nominal_MPs + 6) + 
            log(total_length) + (1 | corral),
          family = poisson(link = "log"),
          ziformula = ~ 1,
          data = GITs)

plot(simulateResiduals(GIT_mod6))  # pass

summary(GIT_mod6)

GIT_mod6.1 <-
  glmmTMB(adjusted_count ~ 1,
          family = poisson(link = "log"),
          ziformula = ~ 1,
          data = GITs)
anova(GIT_mod6, GIT_mod6.1)  # chi2 = 4779, p < 0.001

GIT_mod_pred <- 
  as.data.frame(predict_response(GIT_mod6, 
                                 terms = reference_grid),
                terms_to_colnames = TRUE)
GIT_mod_sim <- 
  as.data.frame(predict_response(GIT_mod6, 
                                 terms = reference_grid,
                                 type = "simulate",
                                 nsim = 1000),
                terms_to_colnames = TRUE)

# Plot model predictions

tiff("GIT Adjusted Data.tiff", width = 18, height = 8, units = "cm",
     res = 500)

set.seed(242)

ggplot(GIT_mod_pred) +
  geom_ribbon(data = GIT_mod_sim,
              aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high),
              alpha = 0.4,
              fill = "grey") +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  geom_point(data = GITs,
             aes(x = nominal_MPs,
                 y = adjusted_count),
             colour = "red",
             size = 1.5,
             alpha = 0.5,
             position = position_jitter(height = 0,
                                        width = 0.1,
                                        seed = 123)) +
  geom_segment(data = GITs,
               aes(x = nominal_MPs,
                   y = count,
                   yend = adjusted_count),
               colour = "red",
               position = position_jitter(height = 0,
                                          width = 0.1,
                                          seed = 123),
               alpha = 0.5,
               linewidth = 0.5) +
  geom_point(data = GITs,
             aes(x = nominal_MPs,
                 y = count,
                 fill = polymer,),
             shape = 21,
             size = 1.5,
             alpha = 0.5,
             position = position_jitter(height = 0,
                                        width = 0.1,
                                        seed = 123)) +
  facet_grid(. ~ polymer,
             scales = "free_y") +
  scale_fill_manual(values = c("pink",
                               "yellow",
                               "blue"),
                    name = "") +
  scale_colour_manual(values = c("pink4",
                                 "yellow4",
                                 "blue4"),
                      name = "") +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs),
                     expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.02, 0),
                     trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000)) +
  labs(x = expression(paste("Nominal Exposure Microplastics "*L^-1)),
       y = expression(paste("Microplastics "*GIT^-1))) +
  theme1 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

###### Inference ######

predict_response(GIT_mod6, terms = reference_grid)

predict_response(GIT_mod6, terms = c("total_length", "nominal_MPs"))

###### livers ####

# Poisson, no interactions

liver_mod1 <-
  glmmTMB(adjusted_count ~ 0 + polymer + log(nominal_MPs + 6) + 
            log(total_length) + (1 | corral) + offset(log(liver_weight)),
          family = poisson(link = "log"),
          data = livers)

plot(simulateResiduals(liver_mod1))  # pass

summary(liver_mod1)

liver_mod1.1 <-
  glmmTMB(adjusted_count ~ offset(log(liver_weight)),
          family = poisson(link = "log"),
          data = livers)

anova(liver_mod1, liver_mod1.1)  # chi2 = 9.20, p = 0.10

liver_mod_pred <- 
  as.data.frame(predict_response(liver_mod1, 
                                 terms = reference_grid,
                                 condition = c(liver_weight = 
                                                 mean(livers$liver_weight))),
                terms_to_colnames = TRUE)

liver_mod_sim <- 
  as.data.frame(predict_response(liver_mod1, 
                                 terms = reference_grid,
                                 type = "simulate",
                                 nsim = 1000,
                                 condition = c(liver_weight =
                                                 mean(livers$liver_weight))),
                terms_to_colnames = TRUE)

# Plot model predictions

tiff("Livers Adjusted Data.tiff", width = 18, height = 8, units = "cm",
     res = 500)

set.seed(242)

ggplot(liver_mod_pred) +
  geom_ribbon(data = liver_mod_sim,
              aes(x = nominal_MPs,
                  ymin = conf.low / mean(livers$liver_weight),
                  ymax = conf.high / mean(livers$liver_weight)),
              alpha = 0.4,
              fill = "grey") +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low / mean(livers$liver_weight),
                  ymax = conf.high / mean(livers$liver_weight),
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted /  mean(livers$liver_weight),
                colour = polymer)) +
  geom_point(data = livers,
             aes(x = nominal_MPs,
                 y = adjusted_count / liver_weight),
             colour = "red",
             size = 1.5,
             alpha = 0.5,
             position = position_jitter(height = 0,
                                        width = 0.1,
                                        seed = 123)) +
  geom_segment(data = livers,
               aes(x = nominal_MPs,
                   y = count / liver_weight,
                   yend = adjusted_count / liver_weight),
               colour = "red",
               position = position_jitter(height = 0,
                                          width = 0.1,
                                          seed = 123),
               alpha = 0.5,
               size = 0.5) +
  geom_point(data = livers,
             aes(x = nominal_MPs,
                 y = count / liver_weight,
                 fill = polymer,),
             shape = 21,
             size = 1.5,
             alpha = 0.5,
             position = position_jitter(height = 0,
                                        width = 0.1,
                                        seed = 123)) +
  facet_grid(. ~ polymer,
             scales = "free_y") +
  scale_fill_manual(values = c("pink",
                               "yellow",
                               "blue"),
                    name = "") +
  scale_colour_manual(values = c("pink4",
                                 "yellow4",
                                 "blue4"),
                      name = "") +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs),
                     expand = c(0.02, 0)) +
  labs(x = expression(paste("Nominal Exposure Microplastics "*L^-1)),
       y = expression(paste("Microplastics "*g^-1))) +
  theme1 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

###### Inference ######

predict_response(liver_mod1, terms = reference_grid,
                 condition = c(liver_weight =
                                 mean(livers$liver_weight)))

###### Muscle ####

# Poisson, no interactions

muscle_mod1 <-
  glmmTMB(adjusted_count ~ 0 + polymer + log(nominal_MPs + 6) + 
            log(total_length) + (1 | corral) + offset(log(fillet_weight)),
          family = poisson(link = "log"),
          data = muscle)

plot(simulateResiduals(muscle_mod1))  # pass

summary(muscle_mod1)

muscle_mod1.1 <-
  glmmTMB(adjusted_count ~ offset(log(fillet_weight)),
          family = poisson(link = "log"),
          data = muscle)

anova(muscle_mod1, muscle_mod1.1)  # chi2 = 9.35, p = 0.10

muscle_mod_pred <- 
  as.data.frame(predict_response(muscle_mod1, 
                                 terms = reference_grid,
                                 condition = c(fillet_weight = 
                                                 mean(muscle$fillet_weight))),
                terms_to_colnames = TRUE)

muscle_mod_sim <- 
  as.data.frame(predict_response(muscle_mod1, 
                                 terms = reference_grid,
                                 type = "simulate",
                                 nsim = 1000,
                                 condition = c(fillet_weight = 
                                                 mean(muscle$fillet_weight))),
                terms_to_colnames = TRUE)

# Plot model predictions

tiff("Muscle Adjusted Data.tiff", width = 18, height = 8, units = "cm",
     res = 500)

set.seed(242)

ggplot(muscle_mod_pred) +
  geom_ribbon(data = muscle_mod_sim,
              aes(x = nominal_MPs,
                  ymin = conf.low / mean(muscle$fillet_weight),
                  ymax = conf.high / mean(muscle$fillet_weight)),
              alpha = 0.4,
              fill = "grey") +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low / mean(muscle$fillet_weight),
                  ymax = conf.high / mean(muscle$fillet_weight),
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted / mean(muscle$fillet_weight),
                colour = polymer)) +
  geom_point(data = muscle,
             aes(x = nominal_MPs,
                 y = adjusted_count / fillet_weight),
             colour = "red",
             size = 1.5,
             alpha = 0.5,
             position = position_jitter(height = 0,
                                        width = 0.1,
                                        seed = 123)) +
  geom_segment(data = muscle,
               aes(x = nominal_MPs,
                   y = count / fillet_weight,
                   yend = adjusted_count / fillet_weight),
               colour = "red",
               position = position_jitter(height = 0,
                                          width = 0.1,
                                          seed = 123),
               alpha = 0.5,
               size = 0.5) +
  geom_point(data = muscle,
             aes(x = nominal_MPs,
                 y = count / fillet_weight,
                 fill = polymer,),
             shape = 21,
             size = 1.5,
             alpha = 0.5,
             position = position_jitter(height = 0,
                                        width = 0.1,
                                        seed = 123)) +
  facet_grid(. ~ polymer,
             scales = "free_y") +
  scale_fill_manual(values = c("pink",
                               "yellow",
                               "blue"),
                    name = "") +
  scale_colour_manual(values = c("pink4",
                                 "yellow4",
                                 "blue4"),
                      name = "") +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs),
                     expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.02, 0),
                     trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000)) +
  labs(x = expression(paste("Nominal Exposure Microplastics "*L^-1)),
       y = expression(paste("Microplastics "*g^-1))) +
  theme1 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

###### Inference ######

predict_response(muscle_mod1, terms = reference_grid,
                 condition = c(fillet_weight =
                                 mean(muscle$fillet_weight)))

#### Explore relationships among organs ####

# Put data into wide form

fish_wide <-
  fish_full_summary2 %>% 
  select(-dummy, -count, -organ_weight, -GIT_weight, -fillet_weight,
         -liver_weight, -gill_weight, -log_organ_weight) %>% 
  pivot_wider(names_from = c(organ),
              values_from = adjusted_count)

##### GIT and liver #####

organmod1 <- glmmTMB(Liver ~ 
                       log(GIT + 1) * polymer +
                       (1 | corral) +
                       log(total_length),
                     family = poisson(link = "log"),
                     data = fish_wide)

plot(simulateResiduals(organmod1, integerResponse = TRUE))

plotResiduals(organmod1, fish_wide$total_length)

summary(organmod1)

organmod2 <- glmmTMB(Liver ~ 
                       (1 | corral),
                     family = poisson(link = "log"),
                     data = fish_wide)

plot(simulateResiduals(organmod2, integerResponse = TRUE))

anova(organmod1, organmod2) # None of the terms had an effect

##### GIT and fillet #####

organmod3 <- glmmTMB(Muscle ~ 
                       log(GIT + 1) * polymer + log(total_length) +
                       (1 | corral),
                     family = poisson(link = "log"),
                     data = fish_wide)

plot(simulateResiduals(organmod3, integerResponse = TRUE))

plotResiduals(organmod2, fish_wide$total_length)

summary(organmod2)

organmod4 <- glmmTMB(Muscle ~ (1 | corral),
                     family = poisson(link = "log"),
                     data = fish_wide)

anova(organmod3, organmod4)  # not sig.

plot(simulateResiduals(organmod4, integerResponse = TRUE))


#### Explore particle size and shape ####

str(fish2)


# Isolate rows with particle measurement

measurements <-
  fish4 %>% 
  filter(!is.na(length) & organ != "Gill" & polymer != "Contamination")

# Statistical differences??

length.mod1 <- lm(log(length) ~ organ, data = measurements)
summary(length.mod1)

plot(simulateResiduals(length.mod1))

length.mod2 <- lm(log(length) ~ 1, data = measurements)

anova(length.mod1, length.mod2)  # not sig.

width.mod1 <- lm(log(width) ~ organ, data = measurements)
summary(width.mod1)

plot(simulateResiduals(width.mod1))

width.mod2 <- lm(log(width) ~ 1, data = measurements)
anova(width.mod1, width.mod2)  # sig.

length.predict <- 
  as.data.frame(predict_response(length.mod1,
                                 terms = "organ"),
                terms_to_colnames = TRUE)

width.predict <-
  as.data.frame(predict_response(width.mod1,
                                 terms = "organ"),
                terms_to_colnames = TRUE)

size.predict <- 
  length.predict %>% 
  mutate(width.predicted = width.predict$predicted,
         width.conf.low = width.predict$conf.low,
         width.conf.high= width.predict$conf.high)

plot(predict_response(width.mod1,
                      terms = "organ"))

# Plot

tiff("Shape Plot.tiff", width = 18, height = 6, units = "cm",
     res = 500)

ggplot(measurements) +
  geom_point(aes(x = width * 1000,
                 y = length * 1000,
                 fill = polymer),
             shape = 21,
             size = 1) +
  geom_point(data = size.predict,
             aes(x = width.predicted * 1000,
                 y = predicted * 1000),
             colour = "red",
             size = 1.5) +
  geom_linerange(data = size.predict,
                 aes(xmin = width.conf.low * 1000,
                     xmax = width.conf.high * 1000,
                     y = predicted * 1000),
                 colour = "red",
                 linewidth = 0.5) +
  geom_linerange(data = size.predict,
                 aes(ymin = conf.low * 1000,
                     ymax = conf.high * 1000,
                     x = width.predicted * 1000),
                 colour = "red",
                 linewidth = 0.5) +
  facet_grid(. ~ organ) +
  labs(x = expression(paste("Width ("*mu*"m)")),
       y = expression(paste("Length ("*mu*"m)"))) +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink"),
                    name = "Polymer") +
  scale_x_continuous(limits = c(0,700),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0,2500),
                     expand = c(0,0)) +
  theme1

dev.off()

