#all plots for paper in one script

#simple SIMM walkthrough--------------
set.seed(123)
y = matrix(c(5, 5.1, 4.7, 3.6, 3.2, 0, -1, -2, -3, -7, 
             3.1, 5.6, 3.6, 4.7, 1.3, 1, -4, -3, -7, -9), 
            ncol = 2)
colnames(y) = c("iso1", "iso2")
mu_s = matrix(c(-10, 0, 10, -10, 10, 0), ncol = 2, nrow = 3)
sigma_s = matrix(c(1, 1, 1, 1, 1, 1), ncol = 2, nrow = 3)
s_names = c("A", "B", "C")
x = c(1.6, 1.7, 2.1, 2.5, 1.1, 3.7, 4.5, 6.8, 7.1, 7.7)

cosimmr_in_1 = cosimmr_load(formula = y ~ x,
                            source_names = s_names,
                            source_means = mu_s,
                            source_sds = sigma_s)


plot(cosimmr_in_1, colour_by_cov = TRUE, cov_name = "x")

cosimmr_out_1 = cosimmr_ffvb(cosimmr_in_1)

summary(cosimmr_out_1, type = "statistics", obs = 1)

plot(cosimmr_out_1, type = c("prop_histogram", "covariates_plot"), obs = 1, cov_name = "x", one_plot = TRUE)

x_pred = data.frame(cov1 = c(3, 5))
pred_out = predict(cosimmr_out_1, x_pred)

summary(pred_out, type = "statistics", obs = c(1,2))


#simulated data Low values--------------

set.seed(123)  # for reproducibility

N <- 50  # number of observations
J <- 2    # number of isotopes
K <- 3    # number of sources
L <- 2    # number of covariates

# Simulating covariates X
X <- matrix(rnorm(N * L), ncol = L)

# Simulating beta coefficients
beta <- matrix(rnorm(K * L, 0, 1), ncol = K)

# Calculating f_ik
f <- X %*% beta

# Softmax transformation for p_ik
p <- exp(f) / rowSums(exp(f))

# Simulating mean_corr and sd_corr
mean_corr <- matrix(runif(K * J, -10, 10), ncol = J, nrow = K)
sd_corr <- matrix(runif(K * J, 0, 2), ncol = J, nrow = K)


# Don't worry about concentration dependence
q <- matrix(1, ncol = J, nrow = K)

# Simulating response y_i
sigma <- runif(J, 0, 0.2)  # Assumed known for simplicity

# Get the overall mean
mu <- y <- sigma_overall <- matrix(0, nrow = N, ncol = J)
for (i in 1:N) {
  for (j in 1:J) {
    mu[i, j] <- sum(p[i, ] * mean_corr[,j])
    sigma_overall[i,j] <-sum(p[i, ]^2 * sd_corr[,j]) + sigma[j]
    y[i, j] <- rnorm(1, mu[i, j], sigma_overall[i,j])
  }
}


cosimmr_test = cosimmr_load(
  formula = y ~ X - 1,
  source_names = c("a", "b", "c"),
  source_means = mean_corr,
  source_sds = sd_corr,
  scale_x = FALSE
)
plot(cosimmr_test, colour_by_cov = TRUE, cov_name = "X")

cosimmr_test_out = cosimmr_ffvb(cosimmr_test,
                                ffvb_control = list(n_output = 3600, S = 1000, P = 50, beta_1 = 0.75, beta_2 = 0.75, tau
                                                                = 500, eps_0 = 0.0011, t_W = 500)
)

plot(cosimmr_test_out, type = "beta_histogram", cov_name = "X")

p1 = posterior_predictive(cosimmr_test_out)

df = as.data.frame(cosimmr_test_out$output$beta)
dim(df)
p <- list(1:ncol(df))
ggplot(df) + geom_histogram(aes(x = df[,1])) + theme_bw() + geom_vline(xintercept = beta[1,1], col = "red")

p[[1]] = ggplot(df) + geom_histogram(aes(x = df[,1])) + 
            theme_bw() + geom_vline(xintercept = beta[1], col = "red") + 
            xlab("beta[1,1]")

p[[2]] = ggplot(df) + geom_histogram(aes(x = df[,2])) + theme_bw() + geom_vline(xintercept = beta[1,2], col = "red")+ xlab("beta[1,2]")
p[[3]] = ggplot(df) + geom_histogram(aes(x = df[,3])) + theme_bw() + geom_vline(xintercept = beta[1,3], col = "red")+ xlab("beta[1,3]")
p[[4]] = ggplot(df) + geom_histogram(aes(x = df[,4])) + theme_bw() + geom_vline(xintercept = beta[2,1], col = "red")+ xlab("beta[2,1]")
p[[5]] = ggplot(df) + geom_histogram(aes(x = df[,5])) + theme_bw() + geom_vline(xintercept = beta[2,2], col = "red")+ xlab("beta[2,2]")
p[[6]] = ggplot(df) + geom_histogram(aes(x = df[,6])) + theme_bw() + geom_vline(xintercept = beta[2,3], col = "red")+ xlab("beta[2,3]")

library(ggpubr)
ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]])


#Medium values ----------------------

set.seed(123)  # for reproducibility

N <- 200  # number of observations
J <- 3    # number of isotopes
K <- 4    # number of sources
L <- 5    # number of covariates

# Simulating covariates X
X <- matrix(rnorm(N * L), ncol = L)

# Simulating beta coefficients
beta <- matrix(rnorm(K * L, 0, 1), ncol = K)

# Calculating f_ik
f <- X %*% beta

# Softmax transformation for p_ik
p <- exp(f) / rowSums(exp(f))

# Simulating mean_corr and sd_corr
mean_corr <- matrix(runif(K * J, -10, 10), ncol = J, nrow = K)
sd_corr <- matrix(runif(K * J, 0, 2), ncol = J, nrow = K)


# Don't worry about concentration dependence
q <- matrix(1, ncol = J, nrow = K)

# Simulating response y_i
sigma <- runif(J, 0, 0.2)  # Assumed known for simplicity ##Change to 0.2

# Get the overall mean
mu <- y <- sigma_overall <- matrix(0, nrow = N, ncol = J)
for (i in 1:N) {
  for (j in 1:J) {
    mu[i, j] <- sum(p[i, ] * mean_corr[,j])
    sigma_overall[i,j] <-sum(p[i, ]^2 * sd_corr[,j]) + sigma[j]
    y[i, j] <- rnorm(1, mu[i, j], sigma_overall[i,j])
  }
}


cosimmr_test2 = cosimmr_load(
  formula = y ~ X - 1,
  source_names = c("a", "b", "c", "d"),
  source_means = mean_corr,
  source_sds = sd_corr,
  scale_x = FALSE
)
plot(cosimmr_test2, tracers = c(1,2))
plot(cosimmr_test2, tracers = c(1,3))
plot(cosimmr_test2, tracers = c(2,3))
cosimmr_test_out2 = cosimmr_ffvb(cosimmr_test2)

p2 = posterior_predictive(cosimmr_test_out2)

#High values ----------------------

set.seed(123)  # for reproducibility

N <- 500  # number of observations
J <- 4    # number of isotopes
K <- 5    # number of sources
L <- 10    # number of covariates

# Simulating covariates X
X <- matrix(rnorm(N * L), ncol = L)

# Simulating beta coefficients
beta <- matrix(rnorm(K * L, 0, 1), ncol = K)

# Calculating f_ik
f <- X %*% beta

# Softmax transformation for p_ik
p <- exp(f) / rowSums(exp(f))

# Simulating mean_corr and sd_corr
mean_corr <- matrix(runif(K * J, -10, 10), ncol = J, nrow = K)
sd_corr <- matrix(runif(K * J, 0, 2), ncol = J, nrow = K)


# Don't worry about concentration dependence
q <- matrix(1, ncol = J, nrow = K)

# Simulating response y_i
sigma <- runif(J, 0, 0.2)  # Assumed known for simplicity

# Get the overall mean
mu <- y <- sigma_overall <- matrix(0, nrow = N, ncol = J)
for (i in 1:N) {
  for (j in 1:J) {
    mu[i, j] <- sum(p[i, ] * mean_corr[,j])
    sigma_overall[i,j] <-sum(p[i, ]^2 * sd_corr[,j]) + sigma[j]
    y[i, j] <- rnorm(1, mu[i, j], sigma_overall[i,j])
  }
}



cosimmr_test3 = cosimmr_load(
  formula = y ~ X - 1,
  source_names = c("a", "b", "c", "d", "e"),
  source_means = mean_corr,
  source_sds = sd_corr,
  scale_x = FALSE
)

cosimmr_test_out3 = cosimmr_ffvb(cosimmr_test3)

p_3 = posterior_predictive(cosimmr_test_out3)

plot(cosimmr_test_out3, type = "beta_histogram", cov_name = "X")

#geese cosimmr--------------
geese_data = cosimmr::geese_data
Groups = as.factor(geese_data$groups)

cosimmr_geese_in = cosimmr_load(
  formula = geese_data$mixtures ~ Groups -1,
  source_names = geese_data$source_names,
  source_means = geese_data$source_means,
  source_sds = geese_data$source_sds,
  correction_means = geese_data$correction_means,
  correction_sds = geese_data$correction_sds,
  concentration_means = geese_data$concentration_means,
  scale_x = FALSE
)

plot(cosimmr_geese_in, colour_by_cov = TRUE, cov_name = "Groups", 
     xlab = expression(paste(delta^13, "C \u2030", sep = "")),
     ylab = expression(paste(delta^15, "N \u2030", sep = "")),
     title = "Isospace plot of Inger et al Geese data")

cosimmr_geese_out = cosimmr_ffvb(cosimmr_geese_in)

plot(cosimmr_geese_out, type = c("prop_histogram", "covariates_plot"), cov_name = "Groups", obs = 1)

geese_post_pred = posterior_predictive(cosimmr_geese_out)

#geese mixsiar--------

library(MixSIAR)

# Load mix data
mix.filename <- system.file("extdata", "geese_consumer.csv", package = "MixSIAR")
mix <- load_mix_data(filename=mix.filename,
                     iso_names=c("d13C","d15N"),
                     factors="Group",
                     fac_random=FALSE,
                     fac_nested=FALSE,
                     cont_effects=NULL)

# Load source data
source.filename <- system.file("extdata", "geese_sources.csv", package = "MixSIAR")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL,
                           conc_dep=TRUE,
                           data_type="means",
                           mix)

# Load discrimination/TDF data
discr.filename <- system.file("extdata", "geese_discrimination.csv", package = "MixSIAR")
discr <- load_discr_data(filename=discr.filename, mix)


# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- FALSE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)


#jags.1 <- run_model(run="test", mix, source, discr, model_filename, alpha.prior=1)
jags.1 <- run_model(run="normal", mix, source, discr, model_filename, alpha.prior=1)

#Plot MixSIAR proportions
source("plot_mixsiar.R")
library(viridis)
plot_mixsiar(jags.1$BUGSoutput$sims.list$p.fac1, source$source_names)

#geese timing comparison --------
library(microbenchmark)
geese_timing = microbenchmark(cosimmr_geese_out = cosimmr_ffvb(cosimmr_geese_in),
                              jags.1 <- run_model(run="normal", mix, source, discr, model_filename, alpha.prior=1),
                              times = 50L)

#isopod cosimmr--------------
isopod_data = cosimmr::iso_data
Site = as.factor(isopod_data$site)
iso_cosimmr = cosimmr_load(
  formula = as.matrix(isopod_data$mixtures) ~ Site -1,
  isopod_data$source_names,
  as.matrix(isopod_data$source_means),
  as.matrix(isopod_data$source_sds),
  as.matrix(isopod_data$correction_means),
  as.matrix(isopod_data$correction_sds),
  scale_x = FALSE
)

plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(1,2))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(1,3))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(1,4))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(1,5))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(1,6))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(1,7))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(1,8))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(2,3))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(2,4))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(2,5))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(2,6))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(2,7))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(3,8))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(3,4))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(3,5))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(3,6))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(3,7))
plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(3,8))


plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(3,5),
     xlab = expression(paste(C[18], "3", omega, "3", sep = "")), # c18.3w3
     ylab = expression(paste(C[20], "4", omega, "6", sep = "")), # c20.4w6
     title = "Isospace plot of Galloway et al Isopod data")

iso_out = cosimmr_ffvb(iso_cosimmr)

plot(iso_out, type = c("prop_histogram", "covariates_plot"), cov_name = "Site")

plot(iso_out, type = "covariates_plot", cov_name = "Site")

iso_post_pred = posterior_predictive(iso_out)

save(iso_out, file = "iso_out.rda")

#isopod mixsiar---------------
# Load mix data
mix.filename <- system.file("extdata", "isopod_consumer.csv", package = "MixSIAR")
mix <- load_mix_data(filename=mix.filename,
                     iso_names=c("c16.4w3","c18.2w6","c18.3w3","c18.4w3","c20.4w6","c20.5w3","c22.5w3","c22.6w3"),
                     factors="Site",
                     fac_random=TRUE,
                     fac_nested=FALSE,
                     cont_effects=NULL)

# Load source data
source.filename <- system.file("extdata", "isopod_sources.csv", package = "MixSIAR")
source <- load_source_data(filename=source.filename,
                           source_factors=NULL,
                           conc_dep=FALSE,
                           data_type="means",
                           mix)

# Load discrimination/TDF data
discr.filename <- system.file("extdata", "isopod_discrimination.csv", package = "MixSIAR")
discr <- load_discr_data(filename=discr.filename, mix)


# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- FALSE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the JAGS model ("test" first, then "normal")
#jags.1 <- run_model(run="test", mix, source, discr, model_filename,alpha.prior=1)

## "normal" took my laptop ~60 minutes to run
jags.1 <- run_model(run="normal", mix, source, discr, model_filename,alpha.prior=1)

plot_mixsiar(jags.1$BUGSoutput$sims.list$p.fac1, source$source_names)

#timing comparison----------
iso_time_out = microbenchmark(jags.1 <- run_model(run="normal", mix, source, discr, model_filename,alpha.prior=1), 
                              iso_out = cosimmr_ffvb(iso_cosimmr), times = 10L)

#Alligator cosimmr--------------

library(cosimmr)
load("alligator_data.rda")


y = alligator_data$mixtures
s_means = alligator_data$source_means
s_sds = alligator_data$source_sds
c_means = alligator_data$TEF_means
c_sds = alligator_data$TEF_sds
s_names = alligator_data$source_names
habitat = alligator_data$habitat
sex = alligator_data$sex
sclass = alligator_data$sclass
Length = alligator_data$length
sex_sclass = alligator_data$sex_sclass

ali_cosimmrmod5 = cosimmr_load(
  formula = as.matrix(y) ~  Length,
  (s_names),
  as.matrix(s_means),
  as.matrix(s_sds),
  as.matrix(c_means),
  as.matrix(c_sds),
  scale_x = TRUE
)

plot(ali_cosimmrmod5, colour_by_cov = TRUE, cov_name = "Length",
     xlab = expression(paste(delta^13, "C \u2030", sep = "")),
ylab = expression(paste(delta^15, "N \u2030", sep = "")),
title = "Isospace plot of Nifong et al Alligator data")

ali_outmod5 = cosimmr_ffvb(ali_cosimmrmod5)

plot(ali_outmod5, type = c("prop_histogram", "covariates_plot"), cov_name = "length")


post_pred_ali = posterior_predictive(ali_outmod5)

plot(geese_in, colour_by_cov = TRUE, cov_name = "Groups",
     xlab = expression(paste(delta^13, "C (‰)", sep = "")),
     ylab = expression(paste(delta^15, "N (‰)", sep = "")),
     title = "Isospace plot of Inger et al Geese data")

plot(iso_cosimmr, colour_by_cov = TRUE, cov_name = "Site", tracers = c(3,5),
     xlab = expression(paste(C[18], "3", omega, "3", sep = "")),
     ylab = expression(paste(C[20], "4", omega, "6", sep = "")),
     title = "Isospace plot of Galloway et al Isopod data")


#Alligator Mixsiar-----------
mix.filename <- system.file("extdata", "alligator_consumer.csv", package = "MixSIAR")
source.filename <- system.file("extdata", "alligator_sources_simplemean.csv", package = "MixSIAR")
discr.filename <- system.file("extdata", "alligator_TEF.csv", package = "MixSIAR")


n.mod <- 8
mix <- vector("list", n.mod) 
mix[[1]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors=NULL,
                          fac_random=NULL,
                          fac_nested=NULL,
                          cont_effects=NULL)
mix[[2]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors="habitat",
                          fac_random=FALSE,
                          fac_nested=FALSE,
                          cont_effects=NULL)
mix[[3]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors="sex",
                          fac_random=FALSE,
                          fac_nested=FALSE,
                          cont_effects=NULL)
mix[[4]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors="sclass",
                          fac_random=FALSE,
                          fac_nested=FALSE,
                          cont_effects=NULL)
mix[[5]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors=NULL,
                          fac_random=NULL,
                          fac_nested=NULL,
                          cont_effects="Length")
mix[[6]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors=c("sex","sclass"),
                          fac_random=c(FALSE,FALSE),
                          fac_nested=c(FALSE,FALSE),
                          cont_effects=NULL)
mix[[7]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors="sex",
                          fac_random=FALSE,
                          fac_nested=FALSE,
                          cont_effects="Length")
mix[[8]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C","d15N"),
                          factors="sex_sclass",
                          fac_random=FALSE,
                          fac_nested=FALSE,
                          cont_effects=NULL)

# Run the models
source <- vector("list", n.mod)
discr <- vector("list", n.mod)
jags.mod <- vector("list", n.mod)

mod = 5 # We only want length

# load source data
source[[mod]] <- load_source_data(filename=source.filename,
                                  source_factors=NULL,
                                  conc_dep=FALSE,
                                  data_type="means",
                                  mix[[mod]])

# load TEF data
discr[[mod]] <- load_discr_data(filename=discr.filename, mix[[mod]])


# Define model structure and write JAGS model file
model_filename <- paste0("MixSIAR_model_", mod, ".txt")
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix[[mod]], source[[mod]])


jags.mod[[mod]] <- run_model(run="short", mix[[mod]], source[[mod]], discr[[mod]], model_filename, alpha.prior=1)


plot_mixsiar(jags.mod[[5]]$BUGSoutput$sims.list$p.ind, source[[5]]$source_names)

#timing comparison-------
output_timing_aligator = microbenchmark(cosimmr = cosimmr_ffvb(ali_cosimmrmod5),
                                        JAGS = run_model(run="short", mix[[5]], source[[5]], discr[[5]], model_filename, alpha.prior=1),
                                        times = 10L)



