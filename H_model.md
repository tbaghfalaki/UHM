---
title: "Hurdle Models Overview"
output: html_document
---

# Hurdle Models Overview

## 1. Hurdle Poisson Model
The model is expressed as:

$$
Y_i \sim \text{ZIP}(\pi_i, \lambda_i)
$$

The probability mass function is given by:

$$
P_{ZT}(Y_i = y | \lambda_i) = \frac{e^{-\lambda_i} \lambda_i^y}{y! \left( 1 - e^{-\lambda_i} \right)},
$$

with the log of the mean parameter defined as:

$$
\log(\lambda_i) = \bm{x}_{i} \bm{\beta}.
$$

---

## 2. Hurdle Negative Binomial Model
This model can be expressed as:

$$
Y_i \sim \text{ZINB}(\pi_i, \lambda_i, r)
$$

The probability mass function is:

$$
P_{ZT}(Y_i = y | \kappa_i, r) = \frac{\Gamma(y + r)}{y! \Gamma(r) (1 - \kappa_i^r)} \kappa_i^r (1 - \kappa_i)^y, \quad y = 1, 2, \ldots
$$

where:

$$
\lambda_i = r \frac{(1 - \kappa_i)}{\kappa_i}.
$$

---

## 3. Hurdle Generalized Poisson Model
This model is defined as:

$$
Y_i \sim \text{ZIGP}(\pi_i, \lambda_i, \phi)
$$

The probability mass function is:

$$
P_{ZT}(Y_i = y |
