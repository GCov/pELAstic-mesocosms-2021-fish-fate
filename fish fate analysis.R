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

# First separate out just the particles

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

rf2  # OOB error rate = 9.64%%

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
                    ntree = 501)

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
  theme_bw()

ggplot(fish3 %>% filter(!is.na(polymer)),
       aes(x = organ,
           fill = count_ID)) +
  geom_bar() +
  facet_wrap(~ fish_ID, scales = "free_y") +
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

fish4_particle <-
  fish4 %>% 
  filter(!is.na(particle))

# Need to make a new data set with every polymer for every sample then
# populate it

fish_new <-
  expand.grid(fish_ID = unique(fish4$fish_ID),
              organ = unique(fish4$organ),
              polymer = unique(fish4$polymer)[-c(1,4)]) %>%
  left_join(fish4 %>% 
              select(fish_ID,
                     organ,
                     sample_ID, 
                     corral, 
                     nominal_MPs, 
                     total_length, 
                     weight, 
                     GIT_weight, 
                     fillet_weight, 
                     liver_weight, 
                     gill_weight) %>% 
              group_by(fish_ID,
                       organ,
                       sample_ID, 
                       corral, 
                       nominal_MPs, 
                       total_length, 
                       weight, 
                       GIT_weight, 
                       fillet_weight, 
                       liver_weight, 
                       gill_weight) %>% 
              summarize(dummy = length(fish_ID)) %>% 
              ungroup(),
            by = c("fish_ID", "organ")) %>% 
  full_join(fish4_particle) %>% 
  filter(polymer != "Contamination")

fish_new$particle[is.na(fish_new$particle)] <- 0

# Create a summary dataset

fish_summary <-
  fish_new %>% 
  group_by(fish_ID,
           organ,
           corral,
           nominal_MPs,
           total_length,
           weight,
           GIT_weight,
           fillet_weight,
           liver_weight,
           gill_weight,
           polymer) %>% 
  summarize(count = sum(particle)) %>% 
  ungroup()

fish_summary$polymer <- as.character(fish_summary$polymer)
fish_summary$polymer <- as.factor(fish_summary$polymer)

str(fish_summary)

# Plot

ggplot(fish_summary,
       aes(x = nominal_MPs,
           y = count,
           fill = polymer)) +
  geom_point(shape = 21) +
  scale_x_continuous(trans = "log1p",
                     breaks = c(0, 1, 100, 1000, 10000)) +
  facet_grid(polymer ~ organ) +
  labs(y = "MPs per ind",
       x = "Nominal Exposure Concentration MPs") +
  scale_fill_manual(values = c("yellow", "blue", "pink"),
                    name = "Polymer") +
  theme_bw()

# Summary stats

fish_totals <-
  fish_summary %>% 
  group_by(fish_ID, organ, nominal_MPs) %>% 
  summarize(total_count = sum(count)) %>% 
  ungroup()

fish_means <-
  fish_totals %>% 
  group_by(organ, nominal_MPs) %>% 
  summarize(mean = mean(total_count),
            sd = sd(total_count)) %>% 
  ungroup()

fish_means
