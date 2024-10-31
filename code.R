# Clear the workspace by removing all objects
rm(list=ls())

# Load necessary libraries
library(UHM)  # For the ZIHR function
library(cplm)  # For models related to claims data

# Set the working directory where the data file is located
setwd("/Users/taban/Desktop/uhm_code/")

# Read the count data from a text file
count <- read.table("cells.txt", header=TRUE)

# Attach the data frame for easier variable access
attach(count)

# Display the names of the columns in the data frame
names(count)

# Create a frequency table for the cells variable
table(cells)

# Calculate and display the mean of 'cells' for different groups
mean_by_smoker <- tapply(cells, smoker, mean)  # Mean by smoking status
mean_by_weight <- tapply(cells, weight, mean)  # Mean by weight
mean_by_sex <- tapply(cells, sex, mean)        # Mean by sex
mean_by_age <- tapply(cells, age, mean)        # Mean by age

# ============================== First Model: Poisson ==============================
# Create a one-hot encoding for the weight variable
weight1 <- nnet::class.ind(weight)

# Combine the original data with the new weight indicators into a new data frame
Data1 <- data.frame(count, weight1)

# Define the response variable and model formula for the count model
modelY <- cells ~ smoker + normal + smoker * normal
n <- length(cells)  # Number of observations

# Define the model for zero inflation
modelZ <- z ~ smoker + normal

# Fit a Zero-Inflated Poisson model
C <- ZIHR(modelY, modelZ,
          data = Data1, 
          n.chains = 2, 
          n.iter = 10000,
          n.burnin = 5000, 
          n.thin = 1, 
          family = "Poisson")

# Summary of the model results
summary_C <- SummaryZIHR(C)

# Print the DIC and LPML values for the model
DIC_C <- C$DIC
LPML_C <- C$LPML

# ============================= Second Model: Negative Binomial =============================
# Fit a Zero-Inflated Negative Binomial model
CNB <- ZIHR(modelY, modelZ,
            data = Data1, 
            n.chains = 2, 
            n.iter = 10000,
            n.burnin = 5000, 
            n.thin = 1, 
            family = "NB")

# Summary of the Negative Binomial model results
summary_CNB <- SummaryZIHR(CNB)

# Print the DIC and LPML values for the Negative Binomial model
DIC_CNB <- CNB$DIC
LPML_CNB <- CNB$LPML

# ============================= Third Model: Logarithmic =============================
# Fit a Zero-Inflated Logarithmic model
CLog <- ZIHR(modelY, modelZ,
             data = Data1, 
             n.chains = 2, 
             n.iter = 10000,
             n.burnin = 5000, 
             n.thin = 1, 
             family = "Logarithmic")

# Summary of the Logarithmic model results
summary_CLog <- SummaryZIHR(CLog)

# Print the DIC and LPML values for the Logarithmic model
DIC_CLog <- CLog$DIC
LPML_CLog <- CLog$LPML

# ============================= Fourth Model: Generalized Poisson =============================
# Fit a Zero-Inflated Generalized Poisson model
CGP <- ZIHR(modelY, modelZ,
            data = Data1, 
            n.chains = 2, 
            n.iter = 10000,
            n.burnin = 5000, 
            n.thin = 1, 
            family = "GP")

# Summary of the Generalized Poisson model results
summary_CGP <- SummaryZIHR(CGP)

# Print the DIC and LPML values for the Generalized Poisson model
DIC_CGP <- CGP$DIC
LPML_CGP <- CGP$LPML

# ============================= Fifth Model: Bell =============================
# Fit a Zero-Inflated Bell model
CBell <- ZIHR(modelY, modelZ,
              data = Data1, 
              n.chains = 2, 
              n.iter = 10000,
              n.burnin = 5000, 
              n.thin = 1, 
              family = "Bell")

# Summary of the Bell model results
summary_CBell <- SummaryZIHR(CBell)

# Print the DIC and LPML values for the Bell model
DIC_CBell <- CBell$DIC
LPML_CBell <- CBell$LPML

# ============================= Insurance Claims Data Analysis =============================
# Load the claims data and subset for relevant cases
da <- subset(AutoClaim, IN_YY == 1)  # Filter claims data based on a specific condition
da <- transform(da, CLM_AMT5 = CLM_AMT5 / 1000, INCOME = INCOME / 10000)  # Scale variables

# Attach the claims data for easier variable access
attach(da)

# Define the model for the claims amount
modelY <- CLM_AMT5 ~ MVR_PTS + AREA
n <- length(CLM_AMT5)  # Number of observations

# Define the model for zero inflation in claims data
modelZ <- z ~ MVR_PTS + INCOME

# ============================= Sixth Model: Gamma =============================
# Fit a Zero-Inflated Gamma model
CGamma <- ZIHR(modelY, modelZ,
               data = da, 
               n.chains = 2, 
               n.iter = 10000,
               n.burnin = 5000, 
               n.thin = 1, 
               family = "Gamma")

# Summary of the Gamma model results
summary_CGamma <- SummaryZIHR(CGamma)

# Print the DIC and LPML values for the Gamma model
DIC_CGamma <- CGamma$DIC
LPML_CGamma <- CGamma$LPML

# ============================= Bayesian Count Data Model with cplm =============================
# Create a new cplm object and examine its structure
obj <- new("cplm")
names(obj)  # Display the names of slots in the cplm object

# Fit a cpglm model
P1 <- cpglm(CLM_AMT5 ~ MVR_PTS + AREA, data = da)

# Summary of the cpglm model results
summary_P1 <- SummaryZIHR(P1)

# Extract coefficients and other parameters from the fitted model
model_parameters <- c(coef(P1), p = P1$p, phi = P1$phi)

# Set random seed for reproducibility
set.seed(10)

# Fit a Bayesian Count Data model with specific settings
B1 <- bcplm(CLM_AMT5 ~ MVR_PTS + AREA, 
            data = da, 
            n.chains = 2, 
            n.iter = 7000, 
            tune.iter = 4000,
            n.burnin = 2000, 
            n.thin = 5, 
            bound.p = c(1.1, 1.95))

# Summary of the Bayesian model results
summary_B1 <- SummaryZIHR(gelman.diag(B1$sims.list)[[1]][, 1])

# Create diagnostic plots for the Bayesian model
xyplot(B1$sims.list[, c(1:2, 20, 21)], xlab = NULL)  # Scatter plot of selected parameters
densityplot(B1$sims.list[, c(1:2, 20, 21)], ylab = NULL)  # Density plot of selected parameters
