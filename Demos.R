library(parallel)
library(dynr)
library(OpenMx)
library(igraph)
library(qgraph)
Trues = readRDS("./TDrifts.RDS")
dataset = data.frame(readRDS("./Dataset.RDS"))
head(dataset)
dataset$id = dataset$ID
source("https://raw.githubusercontent.com/JPark93/ctgimme/refs/heads/main/ctsgimme-ver-0.0.3.R")
# varnames = vector of variable names in dataset used for analysis
# dataframe = the dataset as a data.frame object or as a matrix
# id = Character string denoting the ID variable; ideally, set to "ID"
# time = Character string denoting the time index
  # If you encounter errors, rename the time-column to "Time"
    # This is a bug I'm working on
  # For each subject, time should begin at t = 0
  # Each successive time-point encapsulates time since 0
# ME.var = Measurement error variances
  # Define a p x p matrix of the expected ME variances
  # Defaults to a diagonal matrix with near-zero variances
  # Denotes a model with almost no measurement errors
# PE.var = Process noise variances
  # Define a p x p matrix of the expected process noise variances
  # Defaults to an Identity matrix
# cores = Number of computational cores to use for parallelization
  # WARNING: Do not use more cores than you have
  # Defaults to 4
# directory = Existing folder to save results to
# ben.hoch = Alpha-correction for multiple testing
  # Defaults to TRUE
# Galpha = p-value for group-level effects; when ben.hoch = TRUE, this is the starting value for the first test
  # All other successive tests use a more stringent p-value
# Ialpha = p-value for individual-level effects
# sig.thrsh - Proportion of sample that must have a statistically significant effect for a path to be added to the "group"-level
  # Classic GIMME uses 75% threshold; testing for continuous-time suggests that this can be lower for CT-VAR

example = ctsgimme(varnames = paste0("y", 1:5), 
                   dataframe = dataset,
                   id = "id", 
                   time = "Time",
                   ME.var = diag(1e-5, 5), 
                   PE.var = diag(1.00, 5),
                   cores = 5, 
                   directory = "./ExampleResults/",
                   ben.hoch = TRUE, 
                   Galpha = 0.05, 
                   Ialpha = 0.05,
                   sig.thrsh = 0.51)

# Compare to True Values
par(mfrow = c(1,1))
qgraph(t(Trues[[20]]), diag = TRUE,
       theme = "colorblind",
       layout = "circle",
       edge.labels = round(t(Trues[[20]]), 2))
