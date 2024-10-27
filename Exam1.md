Getting Started
---------------

```
library(UHM)
```
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

Dynamic prediction
---------------
##### By considering current value for the association 

To compute dynamic predictions (DP), various functions are available. 

The initial method involves computing DP through a first-order approximation utilizing the *DP_CV* function, which requires the following arguments:
-  object an object inheriting from class ZIJMCV
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  s the landmark time for prediction
-  t the window of prediction for prediction
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
-  n.thin integer specifying the thinning of the chains; default is 1.

The second approach entails computing DP via an MCMC approximation, enabling the computation of credible intervals. This is achieved using the *DP_CV_CI* function, which necessitates the following arguments:

-  object an object inheriting from class ZIJMCV
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  s the landmark time for prediction
-  t the window of prediction for prediction
-  mi the number of multiple imputation for Monte-Carlo approximation; default is 10.
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
-  n.thin integer specifying the thinning of the chains; default is 1.

The final method involves generating plots for DP using the *DPplot1* function, which requires the following arguments:

-  object an object inheriting from class ZIJMCV
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  s the landmark time for prediction
-  id_new id number for individual who want to plot his/her DP
-  mi the number of multiple imputation for Monte-Carlo approximation; default is 10.
-  by number: increment of the sequence of DP.
-  Marker_lab the label for the response axis
-  Time_lab the label for the time axis
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.


Example: 
```
 DD <- DP_CV(
    object = Z2, s = 0.5, t = 0.5, n.chains = 1, n.iter = 2000, n.burnin = 1000,
    n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )
```

The outputs of this function are as follows:

```
$DP
     id        est
1     5 0.54272784
2    10 0.36723442
3    15 0.21209392
4    16 0.42057850
5    18 0.58191392
6    19 0.93703289
7    21 0.78147087
8    22 0.86166804
.
.
.
144 488 0.79895668
145 492 0.19337471
146 493 0.08209953
147 494 0.42388086
148 495 0.18696540
149 498 0.74811362
150 500 0.66517540

$s
[1] 0.5

$t
[1] 0.5
```
Computing AUC and BS for the predictions
---------------
For this purpose, we use DPCri package <https://github.com/tbaghfalaki/DPCri>.

Computing the criteria using this package is straightforward, as demonstrated by the following commands:

- s the landmark time for prediction
- t the window of prediction for prediction
- Survt the survival time
- CR the indicator for competing risks or censoring
- P the risk predictions
- cause the main cause for prediction


Consider the following command: 

```
library(survival)
library(DPCri)

Criteria(
  s = 0.5, t = 0.5, Survt = dataSurv_v$survtime,
  CR = dataSurv_v$death, P = DD$DP$est, cause = 1
)
```
with the following outputs:

```
$Cri
          est         sd
AUC 0.7843750 0.05783069
BS  0.1673699 0.02236938
```

Additionally, for the plot, we have the following visualization for subject 167:

```
DPplot1(Z2,
  s = 1.5, id_new = 167, by = 0.1, mi = 5,
  Marker_lab="Marker", Time_lab="Time",
  n.chains = 1, n.iter = 20, n.burnin = 10,
  dataLong = dataLong_v, dataSurv = dataSurv_v
 )
```


##### By considering shared random effects

To compute dynamic predictions (DP), various functions are available. 

The initial method involves computing DP through a first-order approximation utilizing the *DP_SRE* function, which requires the following arguments:
-  object an object inheriting from class VS
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  s the landmark time for prediction
-  t the window of prediction for prediction
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
-  n.thin integer specifying the thinning of the chains; default is 1.

The second approach entails computing DP via an MCMC approximation, enabling the computation of credible intervals. This is achieved using the *DP_SRE_CI* function, which necessitates the following arguments:

-  object an object inheriting from class VS
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  s the landmark time for prediction
-  t the window of prediction for prediction
-  mi the number of multiple imputation for Monte-Carlo approximation; default is 10.
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.
-  n.thin integer specifying the thinning of the chains; default is 1.


The final method involves generating plots for DP using the *DPplot2* function, which requires the following arguments:

-  object an object inheriting from class ZISRE
-  dataLong data set of observed longitudinal variables.
-  dataSurv data set of observed survival variables.
-  s the landmark time for prediction
-  id_new id number for individual who want to plot his/her DP
-  mi the number of multiple imputation for Monte-Carlo approximation; default is 10.
-  by number: increment of the sequence of DP.
-  Marker_lab the label for the response axis
-  Time_lab the label for the time axis
-  n.chains the number of parallel chains for the model; default is 1.
-  n.iter integer specifying the total number of iterations; default is 1000.
-  n.burnin integer specifying how many of n.iter to discard as burn-in ; default is 5000.


Example: 

```
DD <- DP_SRE(Z1,
    s = 0.1, t = 0.5, n.chains = 1, n.iter = 2000, n.burnin = 1000,
    n.thin = 1, dataLong = dataLong_v, dataSurv = dataSurv_v
  )

```
with the following outputs:

```
$DP
     id          est
1     5 0.0108342165
2    10 0.8862672150
3    15 0.0315760312
4    16 0.7846448587
5    18 0.1739335618
6    19 0.0342197666
7    21 0.4936661757
8    22 0.8842215307
.
.
.
143 481 0.1513871921
144 488 0.5682667495
145 492 0.1234556872
146 493 0.3496337934
147 494 0.1883084163
148 495 0.6348319831
149 498 0.6761305704
150 500 0.6103851130

$s
[1] 0.1

$t
[1] 0.5
```



Additionally, for the plot, we have the following visualization for subject 498:

```
DPplot2(Z1,
  s = 0.4, id_new = 498, by = 0.2, mi = 5,
  Marker_lab="Biomarker", Time_lab="Time (week)",
  n.chains = 1, n.iter = 20, n.burnin = 10,
  dataLong = dataLong_v, dataSurv = dataSurv_v
)
```




