Getting Started
---------------

```
library(UHM)
```
- Note: In this section, we analyze the dataset dataD from the package, which consists of discrete zero-inflated data. The analysis can be extended to other distributional assumptions by adjusting the family parameter (as used in R's glm() function).


Loading the data from the package includes response variable and covariates. 

```
data(dataD)
```

```
> head(dataD)
  y x1          x2
1 0  0 -0.99818639
2 0  0 -1.28496198
3 0  0 -1.57777270
4 0  0 -0.42469587
5 0  0 -0.06312155
6 1  0  0.74676177
```

The key items here are *modelY* and *modelZ*, representing, respectively, a formula for the mean of the count response and a formula for the probability of a zero outcome. This argument structure is identical to that used in the "glm" function.

In our example, we consider the following models:

```
modelY <- y~x1 + x2
modelZ <- z~x1
```

This versatile package supports a wide range of distributional assumptions, encompassing Gaussian, Gamma, inverse Gaussian, Weibull, exponential, beta, Poisson, negative binomial, logarithmic, Bell, generalized Poisson, and binomial distributions. Further elaboration on the model's specifics can be found in Ganjali et al. (2024).

Finally, we have to use the ZIHR function with the following arguments:

- modelY: a formula for the mean of the count response. This argument is identical to the one in the "glm" function.
- modelZ: a formula for the probability of zero. This argument is identical to the one in the "glm" function.
- data: data set of observed variables.
- n.chains: the number of parallel chains for the model; default is 1.
- n.iter: integer specifying the total number of iterations; default is 1000.
- n.burnin: integer specifying how many of n.iter to discard as burn-in ; default is 5000.
- n.thin: integer specifying the thinning of the chains; default is 1.
- family: Family objects streamline the specification of model details for functions like glm. They cover various distributions like "Gaussian", "Exponential", "Weibull", "Gamma", "Beta", "inverse.gaussian", "Poisson", "NB", "Logarithmic", "Bell", "GP", and "Binomial". Specifically, "NB" and "GP" are tailored for hurdle negative binomial and hurdle generalized Poisson models, respectively, while the others are utilized for the corresponding models based on their names.

As an example, consider the following command, where this implementation has been performed on the data:

```
D1 <- ZIHR(modelY, modelZ,
           data = dataD, n.chains = 2, n.iter = 1000,
           n.burnin = 500, n.thin = 1, family = "Poisson"
)
```

A part of its output of the function is as follows:

```
$modelY
y ~ x1 + x2

$modelZ
z ~ x1

$family
[1] "Poisson"

$MCMC
$MCMC$beta
             [,1]      [,2]     [,3]
   [1,] -1.880414 0.8814448 2.025756
   [2,] -1.893762 0.7522026 2.032844
   [3,] -1.937188 0.7072405 2.006569
   [4,] -1.979348 0.8485283 2.009982
   [5,] -2.010113 1.1346140 2.003349
   [6,] -2.044763 1.1959385 2.011071
   [7,] -2.034832 1.0024188 1.994125
   [8,] -1.997136 1.2474106 1.996354
   [9,] -2.084062 1.2243783 1.983649
  [10,] -2.018987 1.1655343 1.968824
  [11,] -1.925154 0.8792903 2.071911
  [12,] -1.988864 1.0713805 2.015874
  [13,] -1.999956 0.7349105 2.043753
  [14,] -1.950170 1.1359450 2.028198
  [15,] -1.973761 0.9706845 1.952756
  [16,] -1.987063 1.0912046 2.002442
.
.
.
[499,] 1.553220 -0.546168001
 [500,] 1.452789 -0.361618993
 [ reached getOption("max.print") -- omitted 500 rows ]


$Estimation
$Estimation$Count_model
                   Est        SD       L_CI      U_CI     Rhat
(Intercept) -2.0064805 0.2227296 -2.3981608 -1.534287 0.999330
x1           0.9057071 0.1469032  0.6207732  1.172196 1.012919
x2           2.0795586 0.1289667  1.8100930  2.301690 0.999356

$Estimation$Zero_inflated_model
                   Est        SD      L_CI        U_CI     Rhat
(Intercept)  1.4231425 0.1252582  1.185091  1.67389635 1.002524
x1          -0.5429006 0.2560306 -1.038567 -0.03827461 1.007747


$DIC
[1] 764.4662

$LPML
[1] -382.0315
```

To summarize the ZIHR function, we should use *SummaryZIHR*, which takes an object of class ZIHR as its input.

For example, consider the following command:
```
> SummaryZIHR(D1)
```
The output of this function is as follows: 
```
$Estimation
$Estimation$Count_model
                   Est        SD       L_CI      U_CI     Rhat
(Intercept) -2.0064805 0.2227296 -2.3981608 -1.534287 0.999330
x1           0.9057071 0.1469032  0.6207732  1.172196 1.012919
x2           2.0795586 0.1289667  1.8100930  2.301690 0.999356

$Estimation$Zero_inflated_model
                   Est        SD      L_CI        U_CI     Rhat
(Intercept)  1.4231425 0.1252582  1.185091  1.67389635 1.002524
x1          -0.5429006 0.2560306 -1.038567 -0.03827461 1.007747


$DIC
[1] 764.4662

$LPML
[1] -382.0315
```


