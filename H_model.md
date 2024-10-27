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
\log(\lambda_i) = x_{1i} \beta,
$$

$$
logit(\pi_i) = x_{2i} \alpha.
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
P_{ZT}(Y_i = y | \lambda_i, \phi) = \left( \frac{\lambda_i}{1 + \phi \lambda_i} \right)^y \frac{(1 + \phi y)^{y - 1}}{y! \left( 1 - \exp\left(-\frac{\lambda_i}{1 + \phi \lambda_i}\right)\right)} \exp\left( -\frac{\lambda_i (1 + \phi y)}{1 + \phi \lambda_i} \right),
$$

with:

$$
\log(\lambda_i) = x_{1i} \beta.
$$

---

## 4. Hurdle Binomial Model
The model is expressed as:

$$
Y_i \sim \text{ZIBin}(\pi_i, \lambda_i)
$$

The probability function can be expressed as:

$$
P_{ZT}(Y_i = y | n, \lambda_i) = \binom{n}{y} \frac{\lambda_i^y (1 - \lambda_i)^{n - y}}{\left( 1 - (1 - \lambda_i)^{n} \right)},
$$

with:

$$
\text{logit}(\lambda_i) = x_{1i} \beta.
$$

---

## 5. Hurdle Bell Model
This model takes the form:

$$
Y_i \sim \text{ZIBell}(\pi_i, \eta_i)
$$

The probability mass function is:

$$
P_{ZT}(Y_i = y | \eta_i) = \frac{\eta_i^y e^{-e^{\eta_i} + 1} B_y}{y!}, \quad y = 1, 2, \ldots,
$$

where:

$$
B_y = \frac{1}{e} \sum_{k=0}^{\infty} \frac{k^y}{k!}
$$

and the log of the mean parameter is:

$$
\log(\eta_i) = x_{1i} \beta.
$$

---

## 6. Hurdle Logarithmic Model
The model is expressed as:

$$
Y_i \sim \text{ZILog}(\pi_i, \lambda_i)
$$

The probability mass function is given by:

$$
P_{ZT}(Y_i = y | \lambda_i) = \frac{-1}{\log(1 - \lambda_i)} \frac{\lambda_i^y}{y}, \quad y = 1, 2, \ldots.
$$

---

## 7. Hurdle Normal Model
The model is defined as:

$$
Y_i \sim \text{ZIN}(\pi_i, \mu_i, \sigma^2)
$$

The density function is given by:

$$
f_{ZT}(y | \mu_i, \sigma^2) = \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-\frac{(y - \mu_i)^2}{2 \sigma^2}}, \quad y \in \mathbb{R},
$$

where:

$$
\log(\mu_i) = x_{1i} \beta.
$$

---

## 8. Hurdle Gamma Model
This model can be expressed as:

$$
Y_i \sim \text{ZIG}(\pi_i, \mu_i, \sigma)
$$

The density function is given by:

$$
f_{ZT}(y | \lambda_i, \sigma) = \frac{\lambda_i^{\sigma} y^{\sigma - 1} \exp(-\lambda_i y)}{\Gamma(\sigma)}, \quad y > 0,
$$

with:

$$
\log(\lambda_i) = \frac{\sigma}{\mu_i}.
$$

---

## 9. Hurdle Inverse-Gaussian Model
Let:

$$
Y_i \sim \text{ZIIG}(\pi_i, \mu_i, \sigma)
$$

The density function is given by:

$$
f_{ZT}(y | \mu_i, \sigma) = \left( \frac{\lambda}{2 \pi y^3} \right)^{1/2} \exp\left(-\frac{\lambda (y - \mu_i)^2}{2 \mu_i^2 y}\right), \quad y > 0,
$$

where the linear predictor for $\mu_i$ is given as before.

---

## 10. Hurdle Weibull Model
The model can be expressed as:

$$
Y_i \sim \text{ZIW}(\pi_i, \mu_i, \kappa)
$$

The density function is:

$$
f_{ZT}(y | \mu_i, \kappa) = \kappa \mu_i y^{\kappa - 1} \exp(-\mu_i y^{\kappa}), \quad y > 0,
$$

with:

$$
\log(\mu_i) = - (x_{1i} \beta).
$$

---

## 11. Hurdle Exponential Model
The model is defined as:

$$
Y_i \sim \text{ZIE}(\pi_i, \mu_i)
$$

The density function is:

$$
f_{ZT}(y | \mu_i) = \mu_i \exp(-\mu_i y), \quad y > 0,
$$

where the linear predictor for $\mu_i$ is given as before.

---

## 12. Hurdle Beta Model
The model can be expressed as:

$$
Y_i \sim \text{ZIBeta}(\pi_i, \mu_i, \phi)
$$

The density function is:

$$
f_{ZT}(y | \mu_i) = y^{\phi \mu_i} (1 - y)^{\phi (1 - \mu_i)}, \quad 0 < y < 1,
$$

with the linear predictor for $\mu_i$ as follows:

$$
logit(\mu_i) = x_{2i} \beta.
$$
