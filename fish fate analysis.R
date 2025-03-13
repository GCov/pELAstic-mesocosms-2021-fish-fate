#### Load Packages ####

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(randomForest)
library(caret)

#### Load Data and Prepare ####

# First load in the fish data

fish <- read.csv("fish_data.csv",
                 header = TRUE,
                 stringsAsFactors = TRUE)

str(fish)

# Add treatments

treatments <- 
  data.frame(corral = as.factor(c("B", "C", "D", "E", "F", "G", "H", "I")),
             nominal_MPs = as.numeric(c(0, 414, 29240, 100, 6, 7071, 0, 1710)))

fish2 <- 
  fish %>% 
  left_join(treatments,
            by = "corral")

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
  summarize(mean = mean(count),
            sd = sd(count))

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
  scale_fill_manual(values = c("blue", "yellow", "pink")) +
  theme_bw()

# Calculate means

SaR_means <-
  SaR_summary %>% 
  group_by(polymer) %>% 
  summarize(mean = mean(recovery / 10),
            sd = sd(recovery / 10))

#### Modeling ####

# See how organ concentration relate to nominal exposure concentrations



