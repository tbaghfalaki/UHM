This document presents a detailed analysis of count data using various statistical models, including Zero-Inflated models. The analysis focuses on the relationship between cell counts and several predictors such as smoking status, weight, sex, and age.

## Setup

We will start by clearing the workspace, loading necessary libraries, and reading the data.

```{r setup, include=FALSE}
# Clear the workspace by removing all objects
rm(list=ls())

# Load necessary libraries
library(UHM)  # For Zero-Inflated Hazard Regression (ZIHR)
library(nnet)  # For class.ind function
library(cplm)  # For count data models
```

## Data Preparation

Next, we set the working directory and load the count data from a text file named `cells.txt`.

```{r data-preparation}
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

# Display the means
mean_by_smoker
mean_by_weight
mean_by_sex
mean_by_age
```

## Model 1: Poisson Regression

In this section, we will fit a Zero-Inflated Poisson model to the count data.

### Model Specification

We create a one-hot encoding for the weight variable and combine it with the count data.

```{r poisson-model}
# Create a one-hot encoding for the weight variable
weight1 <- nnet::class.ind(weight)

# Combine the original data with the new weight indicators into a new data frame
Data1 <- data.frame(count, weight1)

# Define the response variable and model formula for the count model
modelY <- cells ~ smoker + normal + smoker * normal
n <- length(cells)  # Number of observations

# Define the model for zero inflation
modelZ <- z ~ smoker + normal
```

### Model Fitting

We fit the Zero-Inflated Poisson model and summarize the results.

```{r fit-poisson-model}
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

# Display the summary statistics
summary_C
DIC_C
LPML_C
```

## Model 2: Negative Binomial Regression

Next, we will fit a Zero-Inflated Negative Binomial model.

```{r negative-binomial-model}
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

# Display the summary statistics
summary_CNB
DIC_CNB
LPML_CNB
```

## Model 3: Logarithmic Regression

This section will cover the fitting of a Zero-Inflated Logarithmic model.

```{r logarithmic-model}
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

# Display the summary statistics
summary_CLog
DIC_CLog
LPML_CLog
```

## Model 4: Generalized Poisson Regression

We will now fit a Zero-Inflated Generalized Poisson model.

```{r generalized-poisson-model}
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

# Display the summary statistics
summary_CGP
DIC_CGP
LPML_CGP
```

## Model 5: Bell Regression

In this section, we will fit a Zero-Inflated Bell model.

```{r bell-model}
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

# Display the summary statistics
summary_CBell
DIC_CBell
LPML_CBell
```

## Insurance Claims Data Analysis

We will now analyze a dataset related to insurance claims, focusing on the relationship between claim amounts and various predictors.

### Data Preparation

We will load and preprocess the claims data.

```{r claims-data-preparation}
# Load the claims data and subset for relevant cases
da <- subset(AutoClaim, IN_YY == 1)  # Use data in the Yip and Yau paper
da <- transform(da, CLM_AMT5 = CLM_AMT5 / 1000, INCOME = INCOME / 10000)  # Scale the variables
attach(da)

# Define the response variable and model formula for the claims model
modelY_claims <- CLM_AMT5 ~ MVR_PTS + AREA
n_claims <- length(CLM_AMT5)

# Define the model for zero inflation
modelZ_claims <- z ~ MVR_PTS + INCOME
```

### Gamma Regression

We will fit a Zero-Inflated Gamma model for insurance claims.

```{r gamma-model}
# Fit a Zero-Inflated Gamma model
CGamma <- ZIHR(modelY_claims, modelZ_claims,
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

# Display the summary statistics
summary_CGamma
DIC_CGamma
LPML_CGamma
```

## Conclusion

In this document, we performed a comprehensive analysis of count data using various statistical models. The Zero-Inflated models, including Poisson, Negative Binomial, Logarithmic, Generalized Poisson, and Bell models, were fitted to explore the effects of various predictors on cell counts. Additionally, we analyzed insurance claims data using Zero-Inflated Gamma models and Bayesian Count models. The results provide insights into the factors influencing count data outcomes and highlight the importance of model selection based on data characteristics.
