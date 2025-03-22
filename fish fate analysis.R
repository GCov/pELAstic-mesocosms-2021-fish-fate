#### Load Packages ####

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(randomForest)
library(caret)
library(rstan)
library(tidybayes)
library(posterior)
library(brms)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


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
  scale_fill_manual(values = c("blue", "yellow", "pink")) +
  theme_bw()

# Calculate means

SaR_means <-
  SaR_summary %>% 
  group_by(polymer) %>% 
  summarize(mean_prob_recovery = mean(recovery / 10),
            sd_prob_recovery = sd(recovery / 10))

#### Modeling ####

# See how organ concentration relate to nominal exposure concentrations

# I'm going to try this 3 ways:
# 1) GLM with no correction
# 2) GLM with recovery addition and blank subtraction
# 3) GLM fitted with stan incorporating recovery and contamination

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

##### First method: GLMM #####

GLM1 <- glmmTMB(count ~ 0 + polymer + organ + log(nominal_MPs + 6) +
                  (1 | corral) + (1 | fish_ID),
                family = poisson(link = "log"),
                data = fish_full_summary)

summary(GLM1)

plot(simulateResiduals(GLM1, integerResponse = TRUE))  # looks good!

reference_grid <-
  expand.grid(nominal_MPs = unique(fish_full_summary$nominal_MPs),
              polymer = unique(fish_full_summary$polymer),
              organ = unique(fish_full_summary$organ))

GLM1_pred <- 
  as.data.frame(predict_response(GLM1, 
                                 terms = reference_grid),
                terms_to_colnames = TRUE)

# Plot model predictions

ggplot(GLM1_pred) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  geom_point(data = fish_full_summary,
             aes(x = nominal_MPs,
                 y = count,
                 fill = polymer),
             shape = 21) +
  facet_wrap(polymer ~ organ,
             scales = "free_y") +
  scale_fill_manual(values = c("pink",
                               "yellow",
                               "blue")) +
  scale_colour_manual(values = c("pink",
                                 "yellow",
                                 "blue")) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs)) +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

# The model seems to be under-fitting at the high and low ends for PE
# However DHARMa suggests the fit is fine

##### Second method: GLMM with addition/subtraction #####

levels(SaR_means$polymer) <- c("PET", "PE", "PS")

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
  geom_abline(aes(intercept = 0, slope = 1))

# Now fit the model

GLM2 <- glmmTMB(adjusted_count ~ 0 + polymer + organ + log(nominal_MPs + 6) +
                  (1 | corral / fish_ID),  # Need to figure out if this is OK
                family = poisson(link = "log"),
                data = fish_full_summary2)

summary(GLM2)

plot(simulateResiduals(GLM2))  # looks good!

GLM2_pred <- 
  as.data.frame(predict_response(GLM2, 
                                 terms = reference_grid),
                terms_to_colnames = TRUE)

# Plot model predictions

ggplot(GLM2_pred) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  geom_point(data = fish_full_summary2,
             aes(x = nominal_MPs,
                 y = adjusted_count,
                 fill = polymer),
             shape = 21) +
  facet_wrap(polymer ~ organ,
             scales = "free_y") +
  scale_fill_manual(values = c("pink",
                               "yellow",
                               "blue")) +
  scale_colour_manual(values = c("pink",
                                 "yellow",
                                 "blue")) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs)) +
  scale_y_continuous(trans = "log1p") +
  labs(x = "Nominal MPs per L",
       y = "Adjusted # MPs") +
  theme_bw()

##### Compare the first two models #####

# Compare inference from the two models

GLM1_pred$model <- "Model 1"
GLM2_pred$model <- "Model 2"

both_models <- rbind(GLM1_pred, GLM2_pred)

ggplot(both_models) +
  geom_errorbar(aes(x = as.factor(nominal_MPs),
                    ymin = conf.low,
                    ymax = conf.high,
                    colour = "polymer"),
                position = position_dodge2()) +
  geom_point(aes(x = as.factor(nominal_MPs),
                 y = predicted,
                 fill = model),
             shape = 21,
             position = position_dodge2(width = 0.9)) +
  scale_y_continuous(trans = "log1p") +
  facet_wrap(polymer ~ organ,
             scales = "free_y") +
  scale_fill_manual(values = c("purple4", "red")) +
  scale_colour_manual(values = c("purple4", "red")) +
  labs(x = "Nominal MPs per L",
       y = "Predicted # MPs") +
  theme_bw()

##### Third method: Bayesian GLMM with built-in uncertainty #####

###### Prepare data ######

fish_full_summary$organ <- as.factor(as.character(fish_full_summary$organ))
fish_full_summary$polymer <- as.factor(as.character(fish_full_summary$polymer))
fish_full_summary$fish_ID <- as.factor(as.character(fish_full_summary$fish_ID))

# Min/max scale nominal_MPs

fish_full_summary <-
  fish_full_summary %>% 
  mutate(scaled_MPs = (nominal_MPs - min(nominal_MPs)) / 
           (max(nominal_MPs) - min(nominal_MPs)))

# Create a reference grid for prediction

fish_grid <-
  with(fish_full_summary,
       expand.grid(organ = unique(organ),
                   polymer = unique(polymer),
                   nominal_MPs = seq(0, 30000, length.out = 100))) %>% 
  mutate(scaled_MPs = (nominal_MPs - min(nominal_MPs)) / 
           (max(nominal_MPs) - min(nominal_MPs)))
  

###### Choose priors ######

# Test priors for appropriateness

GLM3_prior <-
  function() {
    beta_nominal_MPs <- rnorm(1, 1, 1)
    beta_organ <- rnorm(length(unique(fish_full_summary$organ)), 0, 1)
    beta_polymer <- rnorm(length(unique(fish_full_summary$polymer)), 0, 1)
    sigma_corral <- rexp(1, 10)
    u_corral <- rnorm(length(unique(fish_full_summary$corral)), 0, sigma_corral)
    sigma_fish <- rexp(1, 10)
    u_fish <- rnorm(length(unique(fish_full_summary$fish_ID)), 0, sigma_fish)
    
    mu <- beta_nominal_MPs * fish_full_summary$scaled_MPs +
      beta_organ[as.numeric(fish_full_summary$organ)] +
      beta_polymer[as.numeric(fish_full_summary$polymer)] +
      u_corral[as.numeric(fish_full_summary$corral)] +
      u_fish[as.numeric(fish_full_summary$fish_ID)]
    
    rpois(n = length(mu), lambda = exp(mu))
  }

set.seed(234624)

prior_sim <-
  replicate(1000, GLM3_prior())

prior_prediction <-
  data.frame(predicted = apply(prior_sim, 1, median),
             conf.low = apply(prior_sim, 1, quantile, probs = 0.025),
             conf.high = apply(prior_sim, 1, quantile, probs = 0.975)) %>% 
  cbind(fish_full_summary)

# Plot the prior prediction
ggplot(prior_prediction) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  facet_wrap(polymer ~ organ) +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow3",
                                 "blue",
                                 "pink")) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs)) +
  scale_y_continuous(trans = "log1p") +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

# I had to tweak some of the values a bit, but these look pretty good now

###### Fit the model #######

blanks_summary$polymer <- factor(blanks_summary$polymer, 
                                 levels = c("PE", "PET", "PS"))
levels(SaR_summary$polymer) <- c("PET", "PE", "PS")
SaR_summary$polymer <- factor(SaR_summary$polymer, levels = c("PE", "PET", "PS"))

GLM3_stan_data <-
  with(fish_full_summary, 
       list(fish_ID = as.integer(fish_ID), 
            n_fish = length(unique(fish_ID)),
            organ = model.matrix(~ organ - 1),
            n_organ = length(unique(organ)),
            polymer = model.matrix(~ polymer - 1),
            n_polymer = length(unique(polymer)),
            corral = as.integer(corral),
            n_corral = length(unique(corral)),
            log_organ_weight = as.numeric(log_organ_weight),
            nominal_MPs = scaled_MPs,
            count = count,
            n = nrow(fish_full_summary),
            blanks_polymer = model.matrix(~ blanks_summary$polymer - 1),
            blanks_count = blanks_summary$count,
            n_blanks = nrow(blanks_summary),
            recovery_polymer = model.matrix(~ SaR_summary$polymer - 1),
            recovery_count = SaR_summary$recovery,
            n_recovery = nrow(SaR_summary),
            n_grid = nrow(fish_grid),
            organ_grid = model.matrix(~ fish_grid$organ - 1),
            polymer_grid = model.matrix(~ fish_grid$polymer - 1),
            MPs_grid = fish_grid$scaled_MPs))

# Start with just a random effects model without blank/recovery
# correction

###### Attempt 1 - No Interactions #######

GLM3_stan_program <- '
data {
  int<lower=1> n;                               // Number of observations
  int<lower=1> n_fish;                          // Number of fish
  int<lower=1> n_corral;                        // Number of corrals
  array[n] int<lower=1, upper=n_fish> fish_ID;  // Fish identity
  array[n] int<lower=1, upper=n_corral> corral; // Corral identiy
  vector[n] nominal_MPs;                        // Continuous effect of nominal MP concentration
  int<lower=1> n_organ;                         // Number of levels in organ
  matrix[n, n_organ] organ;                     // Organ as a matrix of indicators
  int<lower=1> n_polymer;                       // Number of levels in polymer
  matrix[n, n_polymer] polymer;                 // Polymer as a matrix of indicators
  array[n] int<lower=0> count;                  // Microplastic count
  int<lower=1> n_grid;                          // Number of predictions
  matrix[n_grid, n_organ] organ_grid;           // Grid values for organ
  matrix[n_grid, n_polymer] polymer_grid;       // Grid values for polymer
  vector[n_grid] MPs_grid;                      // Grid values for nominal MPs
}

parameters {
  real beta_nominal_MPs;            // Slope for MP concentration
  vector[n_organ] beta_organ;       // Coefficients for organ
  vector[n_polymer] beta_polymer;   // Coefficients for polymer
  vector[n_corral] u_corral;        // Random effect for corral
  real<lower=0> sigma_corral;       // SD of corral random effect
  vector[n_fish] u_fish;            // Random effect for fish
  real<lower=0> sigma_fish;         // SD of fish random effect
}

model {
  // Priors
  beta_nominal_MPs ~ normal(1, 1);
  beta_organ ~ normal(0, 1);
  beta_polymer ~ normal(0, 1);
  u_corral ~ normal(0, sigma_corral);
  sigma_corral ~ exponential(10);
  u_fish ~ normal(0, sigma_fish);
  sigma_fish ~ exponential(10);

  // Likelihood
  count ~ poisson_log(beta_nominal_MPs * nominal_MPs
              + organ * beta_organ
              + polymer * beta_polymer
              + u_corral[corral]
              + u_fish[fish_ID]);
}
 generated quantities {
  real u_fish_sim = normal_rng(0, sigma_fish);
  real u_corral_sim = normal_rng(0, sigma_corral);
 
  array[n_grid] int count_sim =                         // Simulate over grid 
  poisson_log_rng(beta_nominal_MPs * MPs_grid
  + organ_grid * beta_organ 
  + polymer_grid * beta_polymer +
  u_corral_sim                                        // Include random effects
  + u_fish_sim);
  
  array[n] int count_pred =                             // Posterior predictive
  poisson_log_rng(beta_nominal_MPs * nominal_MPs
  + organ * beta_organ 
  + polymer * beta_polymer);
  }
'

set.seed(5454)

GLM3 <- stan(model_code = GLM3_stan_program, data = GLM3_stan_data,
             chains = 4, iter = 5000, warmup = 1000, thin = 1)

traceplot(GLM3)

GLM3_summary <- summarize_draws(GLM3)  # Looks good!

GLM3_draws <-
  GLM3  %>% 
  gather_draws(beta_nominal_MPs, beta_organ[i], beta_polymer[i], 
               sigma_fish, sigma_corral)

GLM3_draws %>%
  ggplot(aes(y = .variable, x = .value, group = i, fill = i)) +
  stat_halfeye(.width = c(.95, .5), alpha = 0.75)

GLM3_draws %>% 
  group_by(.variable,
           i) %>% 
  summarize(mean = mean(.value),
            median = median(.value),
            upper95 = quantile(.value, probs = 0.975),
            lower95 = quantile(.value, probs = 0.025))

# Diagnose with DHARMa

# Extract posterior predictive samples

GLM3_posterior_predict <- t(extract(GLM3)$count_pred)
GLM3_predicted_means <- 
  apply(GLM3_posterior_predict, 1, mean)  # Average predictions for each observation


# Create a DHARMa object for residual diagnostics
GLM3_simulation_output <- createDHARMa(
  simulatedResponse = GLM3_posterior_predict, # Simulated response from posterior predictions
  observedResponse = GLM3_stan_data$count,       # Observed data
  fittedPredictedResponse = GLM3_predicted_means,  # Mean predictions
  integerResponse = TRUE
)

# Plot residual diagnostics
plot(GLM3_simulation_output)

# Check for specific issues
testDispersion(GLM3_simulation_output)  # Test for overdispersion
testZeroInflation(GLM3_simulation_output)  # Test for zero inflation

# Plot simulations from the model

GLM3_simulated <- t(extract(GLM3)$count_sim)

GLM3_predictions <-
  data.frame(predicted = 
               apply(GLM3_simulated, 1, median),
             conf.low = apply(GLM3_simulated, 1, quantile, probs = 0.025),
             conf.high = apply(GLM3_simulated, 1, quantile, probs = 0.975))

GLM3_predictions <-
  fish_grid %>% 
  cbind(GLM3_predictions)

ggplot(GLM3_predictions) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  geom_point(data = fish_full_summary,
             aes(x = nominal_MPs,
                 y = count,
                 fill = polymer),
             shape = 21) +
  facet_wrap(polymer ~ organ,
             scales = "free_y") +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow3",
                                 "blue",
                                 "pink")) +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

ggplot(GLM3_predictions) +
  geom_errorbar(aes(x = as.factor(nominal_MPs),
                    ymin = conf.low,
                    ymax = conf.high,
                    colour = polymer)) +
  geom_point(aes(x = as.factor(nominal_MPs),
                 y = predicted),
             shape = 21) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000)) +
  facet_wrap(polymer ~ organ) +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow",
                                 "blue",
                                 "pink")) +
  labs(x = "Nominal MPs per L",
       y = "Predicted # MPs") +
  theme_bw()

# Maybe needs an interaction? First try a weaker prior on the MP effect to 
# see what happens

###### Attempt 2 - Weak Prior ######

# Try new prior

GLM3.1_prior <-
  function() {
    beta_nominal_MPs <- rnorm(1, 1, 1)
    beta_organ <- rnorm(length(unique(fish_full_summary$organ)), 0, 1)
    beta_polymer <- rnorm(length(unique(fish_full_summary$polymer)), 0, 1)
    sigma_corral <- rexp(1, 10)
    u_corral <- rnorm(length(unique(fish_full_summary$corral)), 0, sigma_corral)
    sigma_fish <- rexp(1, 10)
    u_fish <- rnorm(length(unique(fish_full_summary$fish_ID)), 0, sigma_fish)
    
    mu <- beta_nominal_MPs * fish_full_summary$scaled_MPs +
      beta_organ[as.numeric(fish_full_summary$organ)] +
      beta_polymer[as.numeric(fish_full_summary$polymer)] +
      u_corral[as.numeric(fish_full_summary$corral)] +
      u_fish[as.numeric(fish_full_summary$fish_ID)]
    
    rpois(n = length(mu), lambda = exp(mu))
  }

set.seed(234624)

prior_sim <-
  replicate(1000, GLM3.1_prior())

prior_prediction <-
  data.frame(predicted = apply(prior_sim, 1, median),
             conf.low = apply(prior_sim, 1, quantile, probs = 0.025),
             conf.high = apply(prior_sim, 1, quantile, probs = 0.975)) %>% 
  cbind(fish_full_summary)

# Plot the prior prediction
ggplot(prior_prediction) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  facet_wrap(polymer ~ organ) +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow",
                                 "blue",
                                 "pink")) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs)) +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

GLM3.1_stan_program <- '
data {
  int<lower=1> n;                               // Number of observations
  int<lower=1> n_fish;                          // Number of fish
  int<lower=1> n_corral;                        // Number of corrals
  array[n] int<lower=1, upper=n_fish> fish_ID;  // Fish identity
  array[n] int<lower=1, upper=n_corral> corral; // Corral identiy
  vector[n] nominal_MPs;                        // Continuous effect of nominal MP concentration
  int<lower=1> n_organ;                         // Number of levels in organ
  matrix[n, n_organ] organ;                     // Organ as a matrix of indicators
  int<lower=1> n_polymer;                       // Number of levels in polymer
  matrix[n, n_polymer] polymer;                 // Polymer as a matrix of indicators
  array[n] int<lower=0> count;                  // Microplastic count
  int<lower=1> n_grid;                          // Number of predictions
  matrix[n_grid, n_organ] organ_grid;           // Grid values for organ
  matrix[n_grid, n_polymer] polymer_grid;       // Grid values for polymer
  vector[n_grid] MPs_grid;                     // Grid values for nominal MPs
}

parameters {
  real beta_nominal_MPs;            // Slope for MP concentration
  vector[n_organ] beta_organ;       // Coefficients for organ
  vector[n_polymer] beta_polymer;   // Coefficients for polymer
  vector[n_corral] u_corral;        // Random effect for corral
  real<lower=0> sigma_corral;       // SD of corral random effect
  vector[n_fish] u_fish;            // Random effect for fish
  real<lower=0> sigma_fish;         // SD of fish random effect
}

model {
  // Priors
  beta_nominal_MPs ~ normal(0, 1);
  beta_organ ~ normal(0, 1);
  beta_polymer ~ normal(0, 1);
  u_corral ~ normal(0, sigma_corral);
  sigma_corral ~ exponential(10);
  u_fish ~ normal(0, sigma_fish);
  sigma_fish ~ exponential(10);

  // Likelihood
  count ~ poisson_log(beta_nominal_MPs * nominal_MPs
             + organ * beta_organ
             + polymer * beta_polymer
             + u_corral[corral]
             + u_fish[fish_ID]);
}
 generated quantities {
  array[n_grid] int count_sim =                         // Simulate over grid 
  poisson_log_rng(beta_nominal_MPs * MPs_grid
  + organ_grid * beta_organ 
  + polymer_grid * beta_polymer);
  
  array[n] int count_pred =                                   // Posterior predictive
  poisson_log_rng(beta_nominal_MPs * nominal_MPs
  + organ * beta_organ 
  + polymer * beta_polymer +
  u_corral[corral] +
  u_fish[fish_ID]);
  }
'

set.seed(5454)

GLM3.1 <- stan(model_code = GLM3.1_stan_program, data = GLM3_stan_data,
             chains = 4, iter = 5000, warmup = 1000, thin = 1)

traceplot(GLM3.1)

GLM3.1_summary <- summarize_draws(GLM3.1)  # Looks good!

GLM3.1_draws <-
  GLM3.1  %>% 
  gather_draws(beta_nominal_MPs, beta_organ[i], beta_polymer[i], 
               sigma_fish, sigma_corral)

GLM3.1_draws %>%
  ggplot(aes(y = .variable, x = .value, group = i, fill = i)) +
  stat_halfeye(.width = c(.95, .5), alpha = 0.75)

GLM3.1_draws %>% 
  group_by(.variable,
           i) %>% 
  summarize(mean = mean(.value),
            median = median(.value),
            upper95 = quantile(.value, probs = 0.975),
            lower95 = quantile(.value, probs = 0.025))

# Diagnose with DHARMa

# Extract posterior predictive samples

GLM3.1_posterior_predict <- t(extract(GLM3.1)$count_pred)
GLM3.1_predicted_means <- 
  apply(GLM3.1_posterior_predict, 1, mean)  # Average predictions for each observation


# Create a DHARMa object for residual diagnostics
GLM3.1_simulation_output <- createDHARMa(
  simulatedResponse = GLM3.1_posterior_predict, # Simulated response from posterior predictions
  observedResponse = GLM3_stan_data$count,       # Observed data
  fittedPredictedResponse = GLM3.1_predicted_means  # Mean predictions
)

# Plot residual diagnostics
plot(GLM3.1_simulation_output)  
# Looks like there's some underprediction of higher values

# Check for specific issues
testDispersion(GLM3.1_simulation_output)  # Test for overdispersion
testZeroInflation(GLM3.1_simulation_output)  # Test for zero inflation

# Plot simulations from the model

GLM3.1_simulated <- t(extract(GLM3.1)$count_sim)

GLM3.1_predictions <-
  data.frame(predicted = 
               apply(GLM3.1_simulated, 1, median),
             conf.low = apply(GLM3.1_simulated, 1, quantile, probs = 0.025),
             conf.high = apply(GLM3.1_simulated, 1, quantile, probs = 0.975))

GLM3.1_predictions <-
  fish_grid %>% 
  cbind(GLM3.1_predictions)

ggplot(GLM3.1_predictions) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  geom_point(data = fish_full_summary,
             aes(x = nominal_MPs,
                 y = count,
                 fill = polymer),
             shape = 21) +
  facet_wrap(polymer ~ organ,
             scales = "free_y") +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow",
                                 "blue",
                                 "pink")) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs)) +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

# Setting a weaker prior does appear to strongly affect inference

###### Attempt 3 - Organ:MP Interaction ######

# Try adding an interaction been organ and MPs

GLM3.2_prior <-
  function() {
    beta_nominal_MPs <- rnorm(1, 0, 10)
    beta_organ <- rnorm(length(unique(fish_full_summary$organ)), 0, 10)
    beta_interaction <- rnorm(length(unique(fish_full_summary$organ)), 0, 10)
    beta_polymer <- rnorm(length(unique(fish_full_summary$polymer)), 0, 10)
    sigma_corral <- rexp(1, 10)
    u_corral <- rnorm(length(unique(fish_full_summary$corral)), 0, sigma_corral)
    sigma_fish <- rexp(1, 10)
    u_fish <- rnorm(length(unique(fish_full_summary$fish_ID)), 0, sigma_fish)
    mu <- beta_nominal_MPs * fish_full_summary$scaled_MPs +
      beta_organ[as.numeric(fish_full_summary$organ)] +
      beta_interaction * fish_full_summary$scaled_MPs +
      beta_polymer[as.numeric(fish_full_summary$polymer)] +
      u_corral[as.numeric(fish_full_summary$corral)] +
      u_fish[as.numeric(fish_full_summary$fish_ID)]
    
    rpois(n = length(mu), lambda = exp(mu))
  }

set.seed(234624)

prior_sim <-
  replicate(1000, GLM3.2_prior())

prior_prediction <-
  data.frame(predicted = apply(prior_sim, 1, median),
             conf.low = apply(prior_sim, 1, quantile, probs = 0.025),
             conf.high = apply(prior_sim, 1, quantile, probs = 0.975)) %>% 
  cbind(fish_full_summary)

# Plot the prior prediction
ggplot(prior_prediction) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  facet_wrap(polymer ~ organ) +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow",
                                 "blue",
                                 "pink")) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs)) +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

GLM3.2_stan_program <- '
data {
  int<lower=1> n;                               // Number of observations
  int<lower=1> n_fish;                          // Number of fish
  int<lower=1> n_corral;                        // Number of corrals
  array[n] int<lower=1, upper=n_fish> fish_ID;  // Fish identity
  array[n] int<lower=1, upper=n_corral> corral; // Corral identiy
  vector[n] nominal_MPs;                        // Continuous effect of nominal MP concentration
  int<lower=1> n_organ;                         // Number of levels in organ
  matrix[n, n_organ] organ;                     // Organ as a matrix of indicators
  int<lower=1> n_polymer;                       // Number of levels in polymer
  matrix[n, n_polymer] polymer;                 // Polymer as a matrix of indicators
  array[n] int<lower=0> count;                  // Microplastic count
  int<lower=1> n_grid;                          // Number of predictions
  matrix[n_grid, n_organ] organ_grid;           // Grid values for organ
  matrix[n_grid, n_polymer] polymer_grid;       // Grid values for polymer
  vector[n_grid] MPs_grid;                     // Grid values for nominal MPs
}

parameters {
  vector[n_organ] beta_nominal_MPs; // Slope for MP concentration, varying by organ
  vector[n_organ] beta_organ;       // Coefficients for organ
  vector[n_polymer] beta_polymer;   // Coefficients for polymer
  vector[n_corral] u_corral;        // Random effect for corral
  real<lower=0> sigma_corral;       // SD of corral random effect
  vector[n_fish] u_fish;            // Random effect for fish
  real<lower=0> sigma_fish;         // SD of fish random effect
}

model {
  // Priors
  beta_nominal_MPs ~ normal(0, 0.1);
  beta_organ ~ normal(0, 0.1);
  beta_polymer ~ normal(0, 0.1);
  u_corral ~ normal(0, sigma_corral);
  sigma_corral ~ exponential(10);
  u_fish ~ normal(0, sigma_fish);
  sigma_fish ~ exponential(10);

  // Likelihood
  count ~ poisson_log(organ * beta_organ
             + (organ * beta_nominal_MPs) .* nominal_MPs
             + polymer * beta_polymer
             + u_corral[corral]
             + u_fish[fish_ID]);
}
 generated quantities {
  real u_fish_sim = normal_rng(0, sigma_fish);
  real u_corral_sim = normal_rng(0, sigma_corral);
  
  array[n_grid] int count_sim =                         // Simulate over grid 
  poisson_log_rng(organ_grid * beta_organ 
  + (organ_grid * beta_nominal_MPs) .* MPs_grid
  + polymer_grid * beta_polymer
  + u_fish_sim
  + u_corral_sim);
  
  array[n] int count_pred =                                   // Posterior predictive
  poisson_log_rng(organ * beta_organ
  + (organ * beta_nominal_MPs) .* nominal_MPs
  + polymer * beta_polymer);
  }
'

set.seed(5454)

GLM3.2 <- stan(model_code = GLM3.2_stan_program, data = GLM3_stan_data,
               chains = 4, iter = 10000, warmup = 3000, thin = 1)

traceplot(GLM3.2)

    GLM3.2_summary <- summarize_draws(GLM3.2)  # Looks good!

GLM3.2_draws <-
  GLM3.2  %>% 
  gather_draws(beta_nominal_MPs[i], beta_organ[i], 
               beta_polymer[i], sigma_fish, sigma_corral)

GLM3.2_draws %>%
  ggplot(aes(y = .variable, x = .value, group = i, fill = i)) +
  stat_halfeye(.width = c(.95, .5), alpha = 0.75)

GLM3.2_draws %>% 
  group_by(.variable,
           i) %>% 
  summarize(mean = mean(.value),
            median = median(.value),
            upper95 = quantile(.value, probs = 0.975),
            lower95 = quantile(.value, probs = 0.025))

# Diagnose with DHARMa

# Extract posterior predictive samples

GLM3.2_posterior_predict <- t(extract(GLM3.2)$count_pred)
GLM3.2_predicted_means <- 
  apply(GLM3.2_posterior_predict, 1, mean)  # Average predictions for each observation


# Create a DHARMa object for residual diagnostics
GLM3.2_simulation_output <- createDHARMa(
  simulatedResponse = GLM3.2_posterior_predict, # Simulated response from posterior predictions
  observedResponse = GLM3_stan_data$count,       # Observed data
  fittedPredictedResponse = GLM3.2_predicted_means  # Mean predictions
)

# Plot residual diagnostics
plot(GLM3.2_simulation_output)  
# Looks like there's some underprediction of higher values

# Check for specific issues
testDispersion(GLM3.2_simulation_output)  # Test for overdispersion
testZeroInflation(GLM3.2_simulation_output)  # Test for zero inflation

# Plot simulations from the model

GLM3.2_simulated <- t(extract(GLM3.2)$count_sim)

GLM3.2_predictions <-
  data.frame(predicted = 
               apply(GLM3.2_simulated, 1, median),
             conf.low = apply(GLM3.2_simulated, 1, quantile, probs = 0.025),
             conf.high = apply(GLM3.2_simulated, 1, quantile, probs = 0.975))

GLM3.2_predictions <-
  fish_grid %>% 
  cbind(GLM3.2_predictions)

ggplot(GLM3.2_predictions) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  geom_point(data = fish_full_summary,
             aes(x = nominal_MPs,
                 y = count,
                 fill = polymer),
             shape = 21) +
  facet_wrap(polymer ~ organ,
             scales = "free_y") +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow3",
                                 "blue",
                                 "pink")) +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

# try adding both organ and polymer interactions

###### Attempt 4 - Both Interactions ######

GLM3.3_prior <-
  function() {
    beta_organ <- rnorm(length(unique(fish_full_summary$organ)), 0, 1)
    beta_interaction_organ <- rnorm(length(unique(fish_full_summary$organ)), 
                                    1, 1)
    beta_polymer <- rnorm(length(unique(fish_full_summary$polymer)), 1, 1)
    beta_interaction_polymer <- rnorm(length(unique(fish_full_summary$polymer)), 
                                      1, 1)
    sigma_corral <- rexp(1, 10)
    u_corral <- rnorm(length(unique(fish_full_summary$corral)), 0, sigma_corral)
    sigma_fish <- rexp(1, 10)
    u_fish <- rnorm(length(unique(fish_full_summary$fish_ID)), 0, sigma_fish)
    mu <- beta_organ[as.numeric(fish_full_summary$organ)] +
      beta_interaction_organ[as.numeric(fish_full_summary$organ)] * 
      fish_full_summary$scaled_MPs +
      beta_polymer[as.numeric(fish_full_summary$polymer)] +
      beta_interaction_polymer[as.numeric(fish_full_summary$polymer)] * 
      fish_full_summary$scaled_MPs +
      u_corral[as.numeric(fish_full_summary$corral)] +
      u_fish[as.numeric(fish_full_summary$fish_ID)]
    
    rpois(n = length(mu), lambda = exp(mu))
  }

set.seed(234624)

prior_sim <-
  replicate(1000, GLM3.3_prior())

prior_prediction <-
  data.frame(predicted = apply(prior_sim, 1, median),
             conf.low = apply(prior_sim, 1, quantile, probs = 0.025),
             conf.high = apply(prior_sim, 1, quantile, probs = 0.975)) %>% 
  cbind(fish_full_summary)

# Plot the prior prediction
ggplot(prior_prediction) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  facet_wrap(polymer ~ organ) +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow",
                                 "blue",
                                 "pink")) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs)) +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()


GLM3.3_stan_program <- '
data {
  int<lower=1> n;                               // Number of observations
  int<lower=1> n_fish;                          // Number of fish
  int<lower=1> n_corral;                        // Number of corrals
  array[n] int<lower=1, upper=n_fish> fish_ID;  // Fish identity
  array[n] int<lower=1, upper=n_corral> corral; // Corral identiy
  vector[n] nominal_MPs;                        // Continuous effect of nominal MP concentration
  int<lower=1> n_organ;                         // Number of levels in organ
  matrix[n, n_organ] organ;                     // Organ as a matrix of indicators
  int<lower=1> n_polymer;                       // Number of levels in polymer
  matrix[n, n_polymer] polymer;                 // Polymer as a matrix of indicators
  array[n] int<lower=0> count;                  // Microplastic count
  int<lower=1> n_grid;                          // Number of predictions
  matrix[n_grid, n_organ] organ_grid;           // Grid values for organ
  matrix[n_grid, n_polymer] polymer_grid;       // Grid values for polymer
  vector[n_grid] MPs_grid;                     // Grid values for nominal MPs
}

parameters {
  vector[n_organ] beta_organ;       // Coefficients for organ
  vector[n_organ] beta_nominal_MPs_organ; // Slope for MP concentration, varying by organ
  vector[n_polymer] beta_polymer;   // Coefficients for polymer
  vector[n_polymer] beta_nominal_MPs_polymer; // Slope for MP concentration, varying by polymer
  vector[n_corral] u_corral;        // Random effect for corral
  real<lower=0> sigma_corral;       // SD of corral random effect
  vector[n_fish] u_fish;            // Random effect for fish
  real<lower=0> sigma_fish;         // SD of fish random effect
}

model {
  // Priors
  beta_organ ~ normal(0, 1);
  beta_nominal_MPs_organ ~ normal(1, 1);
  beta_polymer ~ normal(0, 1);
  beta_nominal_MPs_polymer ~ normal(1, 1);
  u_corral ~ normal(0, sigma_corral);
  sigma_corral ~ exponential(10);
  u_fish ~ normal(0, sigma_fish);
  sigma_fish ~ exponential(10);

  // Likelihood
  count ~ poisson_log(organ * beta_organ
             + (organ * beta_nominal_MPs_organ) .* nominal_MPs
             + polymer * beta_polymer
             + (polymer * beta_nominal_MPs_polymer) .* nominal_MPs
             + u_corral[corral]
             + u_fish[fish_ID]);
}
 generated quantities {
  array[n_grid] int count_sim =                         // Simulate over grid 
  poisson_log_rng(organ_grid * beta_organ               // Ignores RFs
  + (organ_grid * beta_nominal_MPs_organ) .* MPs_grid
  + polymer_grid * beta_polymer
  + (polymer_grid * beta_nominal_MPs_polymer) .* MPs_grid);
  
  array[n] int count_pred =                                   // Posterior predictive
  poisson_log_rng(organ * beta_organ
  + (organ * beta_nominal_MPs_organ) .* nominal_MPs
  + polymer * beta_polymer
  + (polymer * beta_nominal_MPs_polymer) .* nominal_MPs
  + u_corral[corral]
  + u_fish[fish_ID]);
  }
'

set.seed(5454)

GLM3.3 <- stan(model_code = GLM3.3_stan_program, data = GLM3_stan_data,
               chains = 4, iter = 10000, warmup = 2000, thin = 1)

traceplot(GLM3.3)

GLM3.3_summary <- summarize_draws(GLM3.3)

GLM3.3_draws <-
  GLM3.3  %>% 
  gather_draws(beta_nominal_MPs_organ[i], beta_organ[i],
               beta_nominal_MPs_polymer[i], beta_polymer[i], 
               sigma_fish, sigma_corral)

GLM3.3_draws %>%
  ggplot(aes(y = .variable, x = .value, group = i, fill = i)) +
  stat_halfeye(.width = c(.95, .5), alpha = 0.75)

GLM3.3_draws %>% 
  group_by(.variable,
           i) %>% 
  summarize(mean = mean(.value),
            median = median(.value),
            upper95 = quantile(.value, probs = 0.975),
            lower95 = quantile(.value, probs = 0.025))

# Diagnose with DHARMa

# Extract posterior predictive samples

GLM3.3_posterior_predict <- t(extract(GLM3.3)$count_pred)
GLM3.3_predicted_means <- 
  apply(GLM3.3_posterior_predict, 1, mean)  # Average predictions for each observation


# Create a DHARMa object for residual diagnostics
GLM3.3_simulation_output <- createDHARMa(
  simulatedResponse = GLM3.3_posterior_predict, # Simulated response from posterior predictions
  observedResponse = GLM3_stan_data$count,       # Observed data
  fittedPredictedResponse = GLM3.3_predicted_means  # Mean predictions
)

# Plot residual diagnostics
plot(GLM3.3_simulation_output)  
# Looks like there's some underprediction of higher values

# Check for specific issues
testDispersion(GLM3.3_simulation_output)  # Test for overdispersion
testZeroInflation(GLM3.3_simulation_output)  # Test for zero inflation

# Plot simulations from the model

GLM3.3_simulated <- t(extract(GLM3.3)$count_sim)

GLM3.3_predictions <-
  data.frame(predicted = 
               apply(GLM3.3_simulated, 1, median),
             conf.low = apply(GLM3.3_simulated, 1, quantile, probs = 0.025),
             conf.high = apply(GLM3.3_simulated, 1, quantile, probs = 0.975))

GLM3.3_predictions <-
  fish_grid %>% 
  cbind(GLM3.3_predictions)

ggplot(GLM3.3_predictions) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  geom_point(data = fish_full_summary,
             aes(x = nominal_MPs,
                 y = count,
                 fill = polymer),
             shape = 21) +
  facet_wrap(polymer ~ organ,
             scales = "free_y") +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow3",
                                 "blue",
                                 "pink")) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs)) +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

# The model is still under-predicting. What if I try the negative binomial 
# distribution? Drop the interactions though.

###### Attempt 5 - Negative Binomial ######
# Negative binomial

GLM3.4_prior <-
  function() {
    beta_organ <- rnorm(length(unique(fish_full_summary$organ)), 0, 10)
    beta_MPs <- rnorm(1, 0, 1)
    beta_polymer <- rnorm(length(unique(fish_full_summary$polymer)), 0, 1)
    sigma_corral <- rexp(1, 10)
    u_corral <- rnorm(length(unique(fish_full_summary$corral)), 0, sigma_corral)
    sigma_fish <- rexp(1, 10)
    u_fish <- rnorm(length(unique(fish_full_summary$fish_ID)), 0, sigma_fish)
    mu <- beta_organ[as.numeric(fish_full_summary$organ)] +
      beta_MPs * fish_full_summary$scaled_MPs +
      beta_polymer[as.numeric(fish_full_summary$polymer)] +
      fish_full_summary$scaled_MPs +
      u_corral[as.numeric(fish_full_summary$corral)] +
      u_fish[as.numeric(fish_full_summary$fish_ID)]
    
    rnbinom(n = length(mu), mu = exp(mu), size = rgamma(1,1,1))
  }

set.seed(234624)

prior_sim <-
  replicate(1000, GLM3.4_prior())

prior_prediction <-
  data.frame(predicted = apply(prior_sim, 1, median),
             conf.low = apply(prior_sim, 1, quantile, probs = 0.025),
             conf.high = apply(prior_sim, 1, quantile, probs = 0.975)) %>% 
  cbind(fish_full_summary)

# Plot the prior prediction
ggplot(prior_prediction) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  facet_wrap(polymer ~ organ) +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow",
                                 "blue",
                                 "pink")) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs)) +
  scale_y_continuous(trans = 'log1p') +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()


GLM3.4_stan_program <- '
data {
  int<lower=1> n;                               // Number of observations
  int<lower=1> n_fish;                          // Number of fish
  int<lower=1> n_corral;                        // Number of corrals
  array[n] int<lower=1, upper=n_fish> fish_ID;  // Fish identity
  array[n] int<lower=1, upper=n_corral> corral; // Corral identiy
  vector[n] nominal_MPs;                        // Continuous effect of nominal MP concentration
  int<lower=1> n_organ;                         // Number of levels in organ
  matrix[n, n_organ] organ;                     // Organ as a matrix of indicators
  int<lower=1> n_polymer;                       // Number of levels in polymer
  matrix[n, n_polymer] polymer;                 // Polymer as a matrix of indicators
  array[n] int<lower=0> count;                  // Microplastic count
  int<lower=1> n_grid;                          // Number of predictions
  matrix[n_grid, n_organ] organ_grid;           // Grid values for organ
  matrix[n_grid, n_polymer] polymer_grid;       // Grid values for polymer
  vector[n_grid] MPs_grid;                     // Grid values for nominal MPs
}

parameters {
  vector[n_organ] beta_organ;       // Coefficients for organ
  real beta_nominal_MPs;            // Slope for MP concentration
  vector[n_polymer] beta_polymer;   // Coefficients for polymer
  vector[n_corral] u_corral;        // Random effect for corral
  real<lower=0> sigma_corral;       // SD of corral random effect
  vector[n_fish] u_fish;            // Random effect for fish
  real<lower=0> sigma_fish;         // SD of fish random effect
  real<lower=0> phi;                // size parameter for NB
}

model {
  // Priors
  beta_organ ~ normal(0, 1);
  beta_nominal_MPs ~ normal(3, 1);
  beta_polymer ~ normal(0, 1);
  u_corral ~ normal(0, sigma_corral);
  sigma_corral ~ exponential(10);
  u_fish ~ normal(0, sigma_fish);
  sigma_fish ~ exponential(10);
  phi ~ gamma(1, 1);

  // Likelihood
  count ~ neg_binomial_2_log(organ * beta_organ
             + beta_nominal_MPs * nominal_MPs
             + polymer * beta_polymer
             + u_corral[corral]
             + u_fish[fish_ID],
             phi);
}
 generated quantities {
  array[n_grid] int count_sim =                         // Simulate over grid 
  neg_binomial_2_log_rng(organ_grid * beta_organ 
  + beta_nominal_MPs * MPs_grid
  + polymer_grid * beta_polymer,
  phi);
  
  array[n] int count_pred =                                   // Posterior predictive
  neg_binomial_2_log_rng(organ * beta_organ
  + beta_nominal_MPs * nominal_MPs
  + polymer * beta_polymer
  + u_corral[corral]
  + u_fish[fish_ID],
  phi);
  }
'

set.seed(5454)

GLM3.4 <- stan(model_code = GLM3.4_stan_program, data = GLM3_stan_data,
               chains = 4, iter = 10000, warmup = 5000, thin = 1)

traceplot(GLM3.4)

GLM3.4_summary <- summarize_draws(GLM3.4)

GLM3.4_draws <-
  GLM3.4  %>% 
  gather_draws(beta_nominal_MPs, beta_organ[i],
               beta_polymer[i], 
               sigma_fish, sigma_corral)

GLM3.4_draws %>%
  ggplot(aes(y = .variable, x = .value, group = i, fill = i)) +
  stat_halfeye(.width = c(.95, .5), alpha = 0.75)

GLM3.4_draws %>% 
  group_by(.variable,
           i) %>% 
  summarize(mean = mean(.value),
            median = median(.value),
            upper95 = quantile(.value, probs = 0.975),
            lower95 = quantile(.value, probs = 0.025))

# Diagnose with DHARMa

# Extract posterior predictive samples

GLM3.4_posterior_predict <- t(extract(GLM3.4)$count_pred)
GLM3.4_predicted_means <- 
  apply(GLM3.4_posterior_predict, 1, mean)  # Average predictions for each observation


# Create a DHARMa object for residual diagnostics
GLM3.4_simulation_output <- createDHARMa(
  simulatedResponse = GLM3.4_posterior_predict, # Simulated response from posterior predictions
  observedResponse = GLM3_stan_data$count,       # Observed data
  fittedPredictedResponse = GLM3.4_predicted_means  # Mean predictions
)

# Plot residual diagnostics
plot(GLM3.4_simulation_output)  
# Looks like there's some underprediction of higher values

# Check for specific issues
testDispersion(GLM3.4_simulation_output)  # Test for overdispersion
testZeroInflation(GLM3.4_simulation_output)  # Test for zero inflation

# Plot simulations from the model

GLM3.4_simulated <- t(extract(GLM3.4)$count_sim)

GLM3.4_predictions <-
  data.frame(predicted = 
               apply(GLM3.4_simulated, 1, median),
             conf.low = apply(GLM3.4_simulated, 1, quantile, probs = 0.025),
             conf.high = apply(GLM3.4_simulated, 1, quantile, probs = 0.975))

GLM3.4_predictions <-
  fish_grid %>% 
  cbind(GLM3.4_predictions)

ggplot(GLM3.4_predictions) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  geom_point(data = fish_full_summary,
             aes(x = nominal_MPs,
                 y = count,
                 fill = polymer),
             shape = 21) +
  facet_wrap(polymer ~ organ,
             scales = "free_y") +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow3",
                                 "blue",
                                 "pink")) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs)) +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

## DHARMa tests pass but getting stan errors and the fitted model doesn't look
## all that much better

## Maybe there needs to be a quadratic function on MPs or something?

# Compare with brms

brm1 <- brm(count ~ 0 + polymer + organ + scaled_MPs:polymer + 
              scaled_MPs:organ +
              (1 | corral) + (1 | fish_ID),
            family = poisson(link = "log"),
            data = fish_full_summary,
            control = list(adapt_delta = 0.99))

summary(brm1)

brm1_model_check <- createDHARMa(
  simulatedResponse = t(posterior_predict(brm1)),
  observedResponse = fish_full_summary$count,
  fittedPredictedResponse = apply(t(posterior_epred(brm1)), 1, mean),
  integerResponse = TRUE)

plot(brm1_model_check)

brm1_posterior <-
  fish_full_summary %>% 
  add_epred_draws(brm1)

ggplot(brm1_posterior) +
  stat_lineribbon(aes(x = nominal_MPs,
                      y = .epred, fill = polymer), 
                  .width = .95,
                  alpha = 0.5,
                  linewidth = 0.5) +
  geom_point(aes(x = nominal_MPs,
                 y = count,
                 fill = polymer),
             shape = 21) +
  facet_wrap(polymer ~ organ,
             scales = "free_y") +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow3",
                                 "blue",
                                 "pink")) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs)) +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

###### Attempt 6 - Corrected ######

####### Just background contamination #######

# Test priors for appropriateness

GLM3.5_prior <-
  function() {
    beta_nominal_MPs <- rnorm(1, 0, 1)
    beta_organ <- rnorm(3, 0, 1)
    beta_polymer <- rnorm(3, 0, 1)
    beta_polymer_blanks <- rnorm(3, 0, 1)
    sigma_fish <- abs(rcauchy(1, 0, 0.1))
    u_fish <- rnorm(length(unique(fish_full_summary$fish_ID)), 0, sigma_fish)
    lambda_contamination <- 
      exp(beta_polymer_blanks[as.numeric(fish_full_summary$polymer)])
    lambda_environment <- 
      exp(beta_nominal_MPs * fish_full_summary$scaled_MPs +
      beta_organ[as.numeric(fish_full_summary$organ)] +
      beta_polymer[as.numeric(fish_full_summary$polymer)] +
      u_fish[as.numeric(fish_full_summary$fish_ID)])
    lambda_total = lambda_environment + lambda_contamination
    rpois(n = length(lambda_total), lambda = lambda_total)
  }

set.seed(6536)

prior_sim <-
  replicate(1000, GLM3.5_prior())

prior_prediction <-
  data.frame(predicted = apply(prior_sim, 1, median),
             conf.low = apply(prior_sim, 1, quantile, probs = 0.025),
             conf.high = apply(prior_sim, 1, quantile, probs = 0.975)) %>% 
  cbind(fish_full_summary)

# Plot the prior prediction
ggplot(prior_prediction) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  facet_wrap(polymer ~ organ) +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow3",
                                 "blue",
                                 "pink")) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs)) +
  scale_y_continuous(trans = "log1p") +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

GLM3.5_stan_program <- '
data {
  int<lower=1> n;                               // Number of observations
  int<lower=1> n_fish;                          // Number of fish
  array[n] int<lower=1, upper=n_fish> fish_ID;  // Fish identity
  vector[n] nominal_MPs;                        // Continuous effect of nominal MP concentration
  int<lower=1> n_organ;                         // Number of levels in organ
  matrix[n, n_organ] organ;                     // Organ as a matrix of indicators
  int<lower=1> n_polymer;                       // Number of levels in polymer
  matrix[n, n_polymer] polymer;                 // Polymer as a matrix of indicators
  array[n] int<lower=0> count;                  // Microplastic count
  int<lower=1> n_blanks;                        // Number of blank values
  matrix[n_blanks, n_polymer] blanks_polymer;   // Polymer for blanks data
  array[n_blanks] int<lower=0> blanks_count;    // Particle counts for blanks
}

parameters {
  vector[n_organ] beta_nominal_MPs;            // Slope for MP concentration, varying by organ
  vector[n_organ] beta_organ;       // Coefficients for organ
  vector[n_polymer] beta_polymer;   // Coefficients for polymer
  vector[n_fish] u_fish;            // Random effect for fish
  real<lower=0> sigma_fish;         // SD of fish random effect
  vector[n_polymer] beta_polymer_blanks; // Coefficient for blanks polymer
}

model {
  // Priors
  beta_nominal_MPs ~ normal(1, 1);
  beta_organ ~ normal(0, 1);
  beta_polymer ~ normal(0, 1);
  u_fish ~ normal(0, sigma_fish);
  sigma_fish ~ cauchy(0, 0.1);
  beta_polymer_blanks ~ normal(0,1);

  // Likelihood
  vector[n_blanks] lambda_contamination = 
                      exp(blanks_polymer * beta_polymer_blanks);
  blanks_count ~ poisson(lambda_contamination);
  
  // Environmental process
  vector[n] lambda_environment = exp((organ * beta_nominal_MPs) .* nominal_MPs
                                      + organ * beta_organ
                                      + polymer * beta_polymer
                                      + u_fish[fish_ID]);
                                      
  // Contamination process
  vector[n] lambda_total = lambda_environment + 
                            exp(polymer * beta_polymer_blanks);
  
  // Observed counts
  count ~ poisson(lambda_total);
}

generated quantities {
  
  // Posterior predictive samples (no random effects)
  array[n] int count_pred;
  vector[n] lambda_environment = exp((organ * beta_nominal_MPs) .* nominal_MPs
                                      + organ * beta_organ 
                                      + polymer * beta_polymer +
                                      u_fish[fish_ID]);
  vector[n] lambda_total = lambda_environment 
                            + exp(polymer * beta_polymer_blanks);
  
  for (i in 1:n) {
    count_pred[i] = poisson_rng(lambda_total[i]);
  }
  }
'

set.seed(5454)

GLM3.5 <- stan(model_code = GLM3.5_stan_program, data = GLM3_stan_data,
             chains = 4, iter = 4000, warmup = 1000, thin = 1)

traceplot(GLM3.5, pars = c("beta_nominal_MPs",
                           "beta_organ",
                           "beta_polymer",
                           "sigma_fish",
                           "beta_polymer_blanks"))

GLM3.5_summary <- summarize_draws(GLM3.5)

GLM3.5_draws <-
  GLM3.5  %>% 
  gather_draws(beta_nominal_MPs[i], 
               beta_organ[i], 
               beta_polymer[i], 
               sigma_fish,
               beta_polymer_blanks[i])

GLM3.5_draws %>%
  ggplot(aes(y = .variable, x = .value, group = i, fill = i)) +
  stat_halfeye(.width = c(.95, .5), alpha = 0.75)

GLM3.5_draws %>% 
  group_by(.variable,
           i) %>% 
  summarize(mean = mean(.value),
            median = median(.value),
            upper95 = quantile(.value, probs = 0.975),
            lower95 = quantile(.value, probs = 0.025))

# Diagnose with DHARMa

# Extract posterior predictive samples

GLM3.5_posterior_predict <- t(extract(GLM3.5)$count_pred)
GLM3.5_predicted_means <- 
  apply(GLM3.5_posterior_predict, 1, mean)  # Average predictions for each observation


# Create a DHARMa object for residual diagnostics
GLM3.5_simulation_output <- createDHARMa(
  simulatedResponse = GLM3.5_posterior_predict, # Simulated response from posterior predictions
  observedResponse = GLM3_stan_data$count,       # Observed data
  fittedPredictedResponse = GLM3.5_predicted_means,  # Mean predictions
  integerResponse = TRUE
)

# Plot residual diagnostics
plot(GLM3.5_simulation_output)

# Check for specific issues
testDispersion(GLM3.5_simulation_output)  # Test for overdispersion
testZeroInflation(GLM3.5_simulation_output)  # Test for zero inflation

# Plot simulations from the model

GLM3.5_simulated <- t(extract(GLM3.5)$count_pred)

GLM3.5_predictions <-
  data.frame(predicted = 
               apply(GLM3.5_simulated, 1, median),
             conf.low = apply(GLM3.5_simulated, 1, quantile, probs = 0.025),
             conf.high = apply(GLM3.5_simulated, 1, quantile, probs = 0.975))

GLM3.5_predictions <-
  fish_full_summary %>% 
  cbind(GLM3.5_predictions)

ggplot(GLM3.5_predictions) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  geom_point(data = fish_full_summary,
             aes(x = nominal_MPs,
                 y = count,
                 fill = polymer),
             shape = 21) +
  facet_wrap(polymer ~ organ,
             scales = "free_y") +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow3",
                                 "blue",
                                 "pink")) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000)) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(GLM3.5_predictions$nominal_MPs)) +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

####### With recovery #######

# Test priors for appropriateness

GLM3.6_prior <-
  function() {
    beta_nominal_MPs <- rnorm(1, 0, 1)
    beta_organ <- rnorm(3, 0, 1)
    beta_polymer <- rnorm(3, 0, 1)
    beta_polymer_blanks <- rnorm(3, 0, 1)
    sigma_fish <- abs(rcauchy(1, 0, 0.1))
    u_fish <- rnorm(length(unique(fish_full_summary$fish_ID)), 0, sigma_fish)
    lambda_contamination <- 
      exp(beta_polymer_blanks[as.numeric(fish_full_summary$polymer)])
    lambda_environment <- 
      exp(beta_nominal_MPs * fish_full_summary$scaled_MPs +
            beta_organ[as.numeric(fish_full_summary$organ)] +
            beta_polymer[as.numeric(fish_full_summary$polymer)] +
            u_fish[as.numeric(fish_full_summary$fish_ID)])
    lambda_total = lambda_environment + lambda_contamination
    rpois(n = length(lambda_total), lambda = lambda_total)
  }

set.seed(6536)

prior_sim <-
  replicate(1000, GLM3.6_prior())

prior_prediction <-
  data.frame(predicted = apply(prior_sim, 1, median),
             conf.low = apply(prior_sim, 1, quantile, probs = 0.025),
             conf.high = apply(prior_sim, 1, quantile, probs = 0.975)) %>% 
  cbind(fish_full_summary)

# Plot the prior prediction
ggplot(prior_prediction) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  facet_wrap(polymer ~ organ) +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow3",
                                 "blue",
                                 "pink")) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(fish_full_summary$nominal_MPs)) +
  scale_y_continuous(trans = "log1p") +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

GLM3.6_stan_program <- '
data {
  int<lower=1> n;                               // Number of observations
  int<lower=1> n_fish;                          // Number of fish
  array[n] int<lower=1, upper=n_fish> fish_ID;  // Fish identity
  vector[n] nominal_MPs;                        // Continuous effect of nominal MP concentration
  int<lower=1> n_organ;                         // Number of levels in organ
  matrix[n, n_organ] organ;                     // Organ as a matrix of indicators
  int<lower=1> n_polymer;                       // Number of levels in polymer
  matrix[n, n_polymer] polymer;                 // Polymer as a matrix of indicators
  array[n] int<lower=0> count;                  // Microplastic count
  int<lower=1> n_blanks;                        // Number of blank values
  matrix[n_blanks, n_polymer] blanks_polymer;   // Polymer for blanks data
  array[n_blanks] int<lower=0> blanks_count;    // Particle counts for blanks
  int<lower=0> n_recovery;                      // Number of recovery samples
  matrix[n_recovery, n_polymer] recovery_polymer; // Polymer for recovery data
  array[n_recovery] int<lower=0> recovery_count;                  // Recovered particles
}

parameters {
  vector[n_organ] beta_nominal_MPs;            // Slope for MP concentration, varying by organ
  vector[n_organ] beta_organ;       // Coefficients for organ
  vector[n_polymer] beta_polymer;   // Coefficients for polymer
  vector[n_fish] u_fish;            // Random effect for fish
  real<lower=0> sigma_fish;         // SD of fish random effect
  vector[n_polymer] beta_polymer_blanks; // Coefficient for blanks polymer
  vector<lower=0, upper=1>[n_polymer] recovery_prob; // Probability of recovery
}

model {
  // Priors
  beta_nominal_MPs ~ normal(1, 1);
  beta_organ ~ normal(0, 1);
  beta_polymer ~ normal(0, 1);
  u_fish ~ normal(0, sigma_fish);
  sigma_fish ~ cauchy(0, 0.1);
  beta_polymer_blanks ~ normal(0,1);
  recovery_prob ~ beta(2, 2);
  
  // Recovery process
  recovery_count ~ binomial(10, recovery_polymer * recovery_prob);

  // Likelihood
  vector[n_blanks] lambda_contamination = 
                      exp(blanks_polymer * beta_polymer_blanks);
  blanks_count ~ poisson(lambda_contamination);
  
  // Environmental process
  vector[n] lambda_environment = exp((organ * beta_nominal_MPs) .* nominal_MPs
                                      + organ * beta_organ
                                      + polymer * beta_polymer
                                      + u_fish[fish_ID]);
                                      
  // Contamination process
  vector[n] lambda_total = lambda_environment + 
                            exp(polymer * beta_polymer_blanks);
  
  // Observed counts
  count ~ poisson(lambda_total .* (polymer * recovery_prob));
}

generated quantities {
  
  // Posterior predictive samples (no random effects)
  array[n] int count_pred;
  array[n] int count_environment;
  vector[n] lambda_environment = exp((organ * beta_nominal_MPs) .* nominal_MPs
                                      + organ * beta_organ 
                                      + polymer * beta_polymer);
  vector[n] lambda_total = lambda_environment 
                            + exp(polymer * beta_polymer_blanks);
                            
  vector[n] lambda_sample = lambda_total .* (polymer * recovery_prob);
  
  for (i in 1:n) {
    count_pred[i] = poisson_rng(lambda_sample[i]);
    count_environment[i] = poisson_rng(lambda_environment[i]);
  }
  }
'

set.seed(5454)

GLM3.6 <- stan(model_code = GLM3.6_stan_program, data = GLM3_stan_data,
               chains = 4, iter = 4000, warmup = 1000, thin = 1)

traceplot(GLM3.6, pars = c("beta_nominal_MPs",
                           "beta_organ",
                           "beta_polymer",
                           "sigma_fish",
                           "beta_polymer_blanks",
                           "recovery_prob"))

GLM3.6_summary <- summarize_draws(GLM3.6)

GLM3.6_draws <-
  GLM3.6  %>% 
  gather_draws(beta_nominal_MPs[i], 
               beta_organ[i], 
               beta_polymer[i], 
               sigma_fish,
               beta_polymer_blanks[i],
               recovery_prob[i])

GLM3.6_draws %>%
  ggplot(aes(y = .variable, x = .value, group = i, fill = i)) +
  stat_halfeye(.width = c(.95, .5), alpha = 0.75)

GLM3.6_draws %>% 
  group_by(.variable,
           i) %>% 
  summarize(mean = mean(.value),
            median = median(.value),
            upper95 = quantile(.value, probs = 0.975),
            lower95 = quantile(.value, probs = 0.025))

# Diagnose with DHARMa

# Extract posterior predictive samples

GLM3.6_posterior_predict <- t(extract(GLM3.6)$count_pred)
GLM3.6_posterior_mean <- t(extract(GLM3.6)$lambda_sample)
GLM3.6_predicted_means <- 
  apply(GLM3.6_posterior_mean, 1, median)  # Average predictions for each observation


# Create a DHARMa object for residual diagnostics
GLM3.6_simulation_output <- createDHARMa(
  simulatedResponse = GLM3.6_posterior_predict, # Simulated response from posterior predictions
  observedResponse = GLM3_stan_data$count,       # Observed data
  fittedPredictedResponse = GLM3.6_predicted_means,  # Mean predictions
  integerResponse = TRUE
)

# Plot residual diagnostics
plot(GLM3.6_simulation_output)

# Check for specific issues
testDispersion(GLM3.6_simulation_output)  # Test for overdispersion
testZeroInflation(GLM3.6_simulation_output)  # Test for zero inflation

# Plot simulations from the model

GLM3.6_simulated <- t(extract(GLM3.6)$count_pred)

GLM3.6_predictions <-
  data.frame(predicted = 
               apply(GLM3.6_simulated, 1, median),
             conf.low = apply(GLM3.6_simulated, 1, quantile, probs = 0.025),
             conf.high = apply(GLM3.6_simulated, 1, quantile, probs = 0.975))

GLM3.6_predictions <-
  fish_full_summary %>% 
  cbind(GLM3.6_predictions)

ggplot(GLM3.6_predictions) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  geom_point(data = fish_full_summary,
             aes(x = nominal_MPs,
                 y = count,
                 fill = polymer),
             shape = 21) +
  facet_wrap(polymer ~ organ,
             scales = "free_y") +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow3",
                                 "blue",
                                 "pink")) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000)) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(GLM3.6_predictions$nominal_MPs)) +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

# Plot lambda_sample with lambda_environment

GLM3.6_sample <- t(extract(GLM3.6)$lambda_sample)

GLM3.6_total_summary <-
  data.frame(predicted.sample = 
               apply(GLM3.6_sample, 1, median),
             conf.low.sample = apply(GLM3.6_sample, 1, quantile, probs = 0.025),
             conf.high.sample = apply(GLM3.6_sample, 1, quantile, probs = 0.975))

GLM3.6_environment <- t(extract(GLM3.6)$lambda_environment)

GLM3.6_environment_summary <-
  data.frame(predicted.environment = 
               apply(GLM3.6_total, 1, median),
             conf.low.environment = apply(GLM3.6_environment, 1, 
                                          quantile, probs = 0.025),
             conf.high.environment = apply(GLM3.6_environment, 1, 
                                           quantile, probs = 0.975))

GLM3.6_real_summary <-
  fish_full_summary %>% 
  cbind(GLM3.6_total_summary) %>% 
  cbind(GLM3.6_environment_summary)
  

ggplot(GLM3.6_real_summary) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low.sample,
                  ymax = conf.high.sample,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted.sample,
                colour = polymer)) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low.environment,
                  ymax = conf.high.environment),
              alpha = 0.3,
              fill = "brown") +
  geom_line(aes(x = nominal_MPs,
                y = predicted.environment),
            colour = "brown") +
  geom_point(data = fish_full_summary,
             aes(x = nominal_MPs,
                 y = count,
                 fill = polymer),
             shape = 21) +
  facet_wrap(polymer ~ organ,
             scales = "free_y") +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow3",
                                 "blue",
                                 "pink")) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000)) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(GLM3.6_predictions$nominal_MPs)) +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()

# Plot simulated 'real' concentrations

GLM3.6_real<- t(extract(GLM3.6)$count_environment)

GLM3.6_real_summary <-
  data.frame(predicted = 
               apply(GLM3.6_real, 1, median),
             conf.low = apply(GLM3.6_real, 1, 
                                          quantile, probs = 0.025),
             conf.high = apply(GLM3.6_real, 1, 
                                           quantile, probs = 0.975))

GLM3.6_real_summary <-
  fish_full_summary %>% 
  cbind(GLM3.6_real_summary)


ggplot(GLM3.6_real_summary) +
  geom_ribbon(aes(x = nominal_MPs,
                  ymin = conf.low,
                  ymax = conf.high,
                  fill = polymer),
              alpha = 0.3) +
  geom_line(aes(x = nominal_MPs,
                y = predicted,
                colour = polymer)) +
  geom_point(data = fish_full_summary,
             aes(x = nominal_MPs,
                 y = count,
                 fill = polymer),
             shape = 21) +
  facet_wrap(polymer ~ organ,
             scales = "free_y") +
  scale_fill_manual(values = c("yellow",
                               "blue",
                               "pink")) +
  scale_colour_manual(values = c("yellow3",
                                 "blue",
                                 "pink")) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000)) +
  scale_x_continuous(trans = "log1p",
                     breaks = unique(GLM3.6_predictions$nominal_MPs)) +
  labs(x = "Nominal MPs per L",
       y = "# MPs") +
  theme_bw()
