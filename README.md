Question 1. In a randomized, double-blind, parallel-group, multicenter study comparing two oral treatments (denoted A and B) for toe-nail infection (De Backer etal., 1998; also see Lesaffre and Spiessons, 2001), patients were evaluated for the degree of onycholysis (the degree of separation of the nail plate from the nail-bed) at baseline (week 0) and at weeks 4, 8, 12, 24, 36, and 48 thereafter. The onycholysis outcome variable is binary (none or mild versus moderate or severe). The binary outcome was evaluated on 294 patients comprising a total of 1908 measurements. The main objective of the analyses is to compare the effects of oral treatments A and B on changes in the probability of the binary onycholysis outcome over the duration of the study. The raw data are stored in an external file: toenail.dat Each row of the data set contains the following five variables: ID,Y,Treatment,Month,Visit. The binary onycholysis outcome variable Y is coded 0 = none or mild, 1 = moderate or severe. The categorical variable Treatment is coded 1=oral treatment A, 0=oral treatment B. The variable Month denotes the exact timing of measurements in months. The variable Visit denotes the visit number (visit numbers 1-7 correspond to scheduled visits at 0, 4, 8, 12, 24, 36, and 48 weeks).
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### 1. Consider a random effects model with a random intercept for the log odds of moderate or severe onycholysis. Assuming linear trends and month as the time variable.

Let *y*<sub>*i**j*</sub> = the binary outcome of the severity of
onycholysis, then we assume
*y*<sub>*i**j*</sub> ∼ *b**i**n**o**m**i**a**l*(*n*, *p*<sub>*i**j*</sub>).

Variance function:
*V**a**r*(*y*<sub>*i**j*</sub>|*b*<sub>*i*</sub>)=*ϕ**v*(*E*\[*y*<sub>*i**j*</sub>|*b*<sub>*i*</sub>\]) = 1 ⋅ *E*\[*y*<sub>*i**j*</sub>|*b*<sub>*i*</sub>\]⋅(1 − *E*\[*y*<sub>*i**j*</sub>|*b*<sub>*i*</sub>\]).

We removed 5 subjects that has only baseline observation. Then add
baseline as another feature, so the model becomes:

$$
\\begin{split}
g(E\[y\_{ij} | b\_{i}\]) 
&= logit(E\[Y\_{ij}|b\_{i}\]) \\\\
&= (\\beta\_{1} + b\_{i}) + \\beta\_{2}Treatment\_{i} + \\beta\_{3}Month\_{ij} + \\beta\_{4}Baseline\_{i} + \\beta\_{5}Treatment\_{i} \* Month\_{ij} \\\\
\\end{split}
$$

where $b\_{i} \\sim MVN(\\underline{0}, G)$, *G* = (*g*<sub>11</sub>)
and *b*<sub>*i*</sub> is independent of X.

### 2. Provide Interpretations for the fixed effects coefficients in your model. Interpret the random effect parameter.

For fixed effects:

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">Estimate</th>
<th align="right">Std. Error</th>
<th align="right">z value</th>
<th align="right">Pr(&gt;|z|)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>(Intercept)</td>
<td align="right">-3.7364682</td>
<td align="right">0.4720985</td>
<td align="right">-7.9145949</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>Treatment1</td>
<td align="right">-0.2127069</td>
<td align="right">0.4788883</td>
<td align="right">-0.4441681</td>
<td align="right">0.6569211</td>
</tr>
<tr class="odd">
<td>Month</td>
<td align="right">-0.3834852</td>
<td align="right">0.0489840</td>
<td align="right">-7.8287900</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>Baseline1</td>
<td align="right">5.7745242</td>
<td align="right">0.5724672</td>
<td align="right">10.0870823</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>Treatment1:Month</td>
<td align="right">-0.1242952</td>
<td align="right">0.0721727</td>
<td align="right">-1.7221905</td>
<td align="right">0.0850350</td>
</tr>
</tbody>
</table>

-   *β*<sub>1</sub>: the log odds of outcome at month 0 is -3.736, for a
    typical individual who belongs to treatment B and whose baseline
    outcome is 0.

-   *β*<sub>2</sub>: the log odds ratio of outcome at month 0 between
    two individuals, who belong to treatment A and treatment B
    respectively, but have similar underline propensity for response and
    have the same outcome at baseline, is -0.213.

-   *β*<sub>3</sub>: the log odds ratio of outcome for 1 unit increase
    in Month is -0.383, given a specific individual belongs to
    treatment B.

-   *β*<sub>4</sub>: the log odds ratio of outcome at a specific month
    between two individuals, who have different outcomes at baseline,
    but have similar underline propensity for response and belong to the
    same group, is 5.775.

-   *β*<sub>3</sub> + *β*<sub>5</sub>: the log odds ratio of outcome for
    1 unit increase in Month is -0.508, given a specific individual
    belongs to treatment A.

-   *β*<sub>5</sub>: the difference in 'log odds ratio of outcome for 1
    unit increase in Month' between two individuals, who belong to
    treatment A and treatment B respectively, is -0.124

For random effect:

<table>
<thead>
<tr class="header">
<th align="left">grp</th>
<th align="left">var1</th>
<th align="left">var2</th>
<th align="right">vcov</th>
<th align="right">sdcor</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Subject_ID</td>
<td align="left">(Intercept)</td>
<td align="left">NA</td>
<td align="right">6.39713</td>
<td align="right">2.529255</td>
</tr>
</tbody>
</table>

Var(*b*<sub>*i*</sub>) = *g*<sub>11</sub> = 6.397, is the variance of
intercept(*β*<sub>1</sub>) random effect.

### 3. From the results of your analysis what conclusions do you draw about the effect of treatment on changes in the severity of onycholysis over time? Provide results that support your conclusions.

Get the p-value for *β*<sub>2</sub> and *β*<sub>5</sub>:
*p*<sub>*β*<sub>2</sub></sub> = 0.657, *p*<sub>*β*<sub>5</sub></sub> =
0.085. All greater than 0.05, so we accept the null hypothesis and state
that there is enough evidence to support 'Treatment' has no effect on
changes in the severity of onycholysis over time.

### 4. How are the interpretations different from the GEE model.

-   GEE interpret the parameters as population average

-   GLMM interpret the parameters as subject-specific

Question 2. The Skin Cancer Prevention Study was a randomized, double-blind, placebo-controlled clinical trial of beta carotene to prevent non-melanoma skin cancer in high-risk subjects (Greenberg et al., 1989, 1990; also see Stukel, 1993). A total of 1805 subjects were randomized to either placebo or 50 mg of beta carotene per day for 5 years. Subjects were examined once a year and biopsied if a cancer was suspected to determine the number of new skin cancers occurring since the last exam. The outcome variable is a count of the number of new skin cancers per year. The outcome was evaluated on 1683 subjects comprising a total of 7081 measurements. The main objective of the analyses is to compare the effects of beta carotene on skin cancer rates. The raw data are stored in an external file: skin.dat Each row of the data set contains the following 9 variables: ID,Center, Age, Skin, Gender, Exposure, Y, Treatment and Year.
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Note: The outcome variable Y is a count of the of the number of new skin cancers per year. The categorical variable Treatment is coded 1=beta carotene, 0 =placebo. The variable Year denotes the year of follow-up. The categorical variable Gender is coded 1 male, 0 female. The categorical variable Skin denotes skin type and is coded 1 = burns, 0 otherwise. The variable Exposure is a count of the number of previous skin cancers. The variable Age is the age (in years) of each subject at randomization.
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### 1. Set up a suitable random effects (random intercept) model for rate of skin cancers with Treatment and Year as covariates.

Let *y*<sub>*i**j*</sub> = the count of new skin cancer per year, then
we assume *y*<sub>*i**j*</sub> ∼ *P**o**s*(*λ*<sub>*i**j*</sub>).

Variance function:
*V**a**r*(*y*<sub>*i**j*</sub>|*b*<sub>*i*</sub>)=*ϕ**v*(*E*\[*y*<sub>*i**j*</sub>|*b*<sub>*i*</sub>\]) = 1 ⋅ *E*\[*y*<sub>*i**j*</sub>|*b*<sub>*i*</sub>\].

the model looks like:

$$
\\begin{split}
g(E\[y\_{ij} | b\_{i}\]) 
&= log(E\[Y\_{ij}|b\_{i}\]) \\\\
&= (\\beta\_{1} + b\_{i}) + \\beta\_{2}Treatment\_{i} + \\beta\_{3}Year\_{ij} + \\beta\_{4}Treatment\_{i} \* Year\_{ij} \\\\
\\end{split}
$$

where $b\_{i} \\sim MVN(\\underline{0}, G)$, *G* = (*g*<sub>11</sub>)
and *b*<sub>*i*</sub> is independent of X.

### 2. Provide Interpretations for the fixed effects coefficients in your model. Interpret the random effect parameter.

For fixed effects:

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">Estimate</th>
<th align="right">Std. Error</th>
<th align="right">z value</th>
<th align="right">Pr(&gt;|z|)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>(Intercept)</td>
<td align="right">-2.4149720</td>
<td align="right">0.1108829</td>
<td align="right">-21.7794798</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>Treatment1</td>
<td align="right">0.0798434</td>
<td align="right">0.1383237</td>
<td align="right">0.5772218</td>
<td align="right">0.5637897</td>
</tr>
<tr class="odd">
<td>Year</td>
<td align="right">-0.0002495</td>
<td align="right">0.0261529</td>
<td align="right">-0.0095418</td>
<td align="right">0.9923869</td>
</tr>
<tr class="even">
<td>Treatment1:Year</td>
<td align="right">0.0348951</td>
<td align="right">0.0358871</td>
<td align="right">0.9723568</td>
<td align="right">0.3308731</td>
</tr>
</tbody>
</table>

-   *β*<sub>1</sub>: the log rate of outcome at year 0 is -2.415, for a
    typical individual who belongs to placebo group.

-   *β*<sub>2</sub>: the log rate ratio of outcome at year 0 between two
    individuals, who belong to beta carotene and placebo respectively,
    but have similar underline propensity for response, is 0.08.

-   *β*<sub>3</sub>: the log rate ratio of outcome for 1 unit increase
    in year is 0, given a specific individual belongs to placebo group.

-   *β*<sub>3</sub> + *β*<sub>4</sub>: the log rate ratio of outcome for
    1 unit increase in year is 0.035, given a specific individual
    belongs to beta carotene.

-   *β*<sub>4</sub>: the difference in 'log rate ratio of outcome for 1
    unit increase in year' between two individuals, who belong to beta
    carotene and placebo respectively, is 0.111

For random effect:

<table>
<thead>
<tr class="header">
<th align="left">grp</th>
<th align="left">var1</th>
<th align="left">var2</th>
<th align="right">vcov</th>
<th align="right">sdcor</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">ID</td>
<td align="left">(Intercept)</td>
<td align="left">NA</td>
<td align="right">2.189255</td>
<td align="right">1.479613</td>
</tr>
</tbody>
</table>

Var(*b*<sub>*i*</sub>) = *g*<sub>11</sub> = 2.189, is the variance of
intercept(*β*<sub>1</sub>) random effect.

### 3. From the results of your analysis what conclusions do you draw about the effect of beta carotene on the rate of skin cancers? Provide results that support your conclusions.

Get the p-value for *β*<sub>2</sub> and *β*<sub>4</sub>:
*p*<sub>*β*<sub>2</sub></sub> = 0.564, *p*<sub>*β*<sub>4</sub></sub> =
0.331. All greater than 0.05, so we accept the null hypothesis and state
that there is enough evidence to support 'beta carotene' has no effect
on changes of the rate of skin cancers.

### 4. Repeat the above analysis adjusting for skin type, age, and the count of the number of previous skin cancers. What conclusions do you draw about the effect of beta carotene on the adjusted rate of skin cancers?

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">Estimate</th>
<th align="right">Std. Error</th>
<th align="right">z value</th>
<th align="right">Pr(&gt;|z|)</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>(Intercept)</td>
<td align="right">-4.1672145</td>
<td align="right">0.3205286</td>
<td align="right">-13.0010699</td>
<td align="right">0.0000000</td>
</tr>
<tr class="even">
<td>Treatment1</td>
<td align="right">0.0268258</td>
<td align="right">0.1295852</td>
<td align="right">0.2070125</td>
<td align="right">0.8360001</td>
</tr>
<tr class="odd">
<td>Year</td>
<td align="right">0.0008950</td>
<td align="right">0.0260899</td>
<td align="right">0.0343055</td>
<td align="right">0.9726335</td>
</tr>
<tr class="even">
<td>Skin1</td>
<td align="right">0.3286485</td>
<td align="right">0.0878530</td>
<td align="right">3.7408891</td>
<td align="right">0.0001834</td>
</tr>
<tr class="odd">
<td>Age</td>
<td align="right">0.0184366</td>
<td align="right">0.0046315</td>
<td align="right">3.9807193</td>
<td align="right">0.0000687</td>
</tr>
<tr class="even">
<td>Exposure</td>
<td align="right">0.1885961</td>
<td align="right">0.0106753</td>
<td align="right">17.6665325</td>
<td align="right">0.0000000</td>
</tr>
<tr class="odd">
<td>Treatment1:Year</td>
<td align="right">0.0367268</td>
<td align="right">0.0357915</td>
<td align="right">1.0261302</td>
<td align="right">0.3048302</td>
</tr>
</tbody>
</table>

Get the p-value for *β*<sub>2</sub> and *β*<sub>4</sub>:
*p*<sub>*β*<sub>2</sub></sub> = 0.836, *p*<sub>*β*<sub>4</sub></sub> =
0.305. All greater than 0.05, so we accept the null hypothesis and state
that there is enough evidence to support 'beta carotene' has no effect
on changes of the rate of skin cancers.

### 5. How are the interpretations different from the GEE model.

-   GEE interpret the parameters as population average

-   GLMM interpret the parameters as subject-specific

Appendix: code
--------------

    knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, comment = "")
    library(tidyverse)
    library(lme4)  # for glmer
    ## Question 1
    # load original data
    toenail1 <- read.csv("toenail.dat", header = T, sep = ' ') %>% 
      as_tibble() %>% 
      mutate(Response = as.factor(Response),
             Treatment = as.factor(Treatment))

    # remove subject only have 1 observation
    toenail2 <- toenail1 %>% 
      filter(!Subject_ID %in% names(which(table(toenail1$Subject_ID) == 1)))

    # add baseline column
    toenail3 <- toenail2 %>% filter(Month != 0)
    toenail3$Baseline <- rep(subset(toenail2, Month == 0)$Response, as.numeric(table(toenail3$Subject_ID)))
      
    toenail_glmm <- glmer(Response ~ Treatment + Month + Baseline + (1 | Subject_ID) + Treatment * Month, 
                          data = toenail3, family = "binomial")

    # fixed effects
    toenail_sum = summary(toenail_glmm)
    toenail_sum$coefficients %>% knitr::kable()

    # random effect
    toenail_sum$varcor %>% knitr::kable()

    ## Question 2
    # load data
    skin <- read.table("skin2.txt", header = F)
    colnames(skin) = c("ID", "Center", "Age", "Skin", "Gender", "Exposure", "Y", "Treatment", "Year")
    skin = skin %>% 
      as_tibble() %>% 
      mutate(Skin = as.factor(Skin),
             Gender = as.factor(Gender),
             Treatment = as.factor(Treatment))

    skin_glmm <- glmer(Y ~ Treatment*Year + (1 | ID), data = skin, family = "poisson")

    # fixed effects
    skin_sum = summary(skin_glmm)
    skin_sum$coefficients %>% knitr::kable()

    # random effect
    skin_sum$varcor %>% knitr::kable()

    # adjusted model
    skin2_glmm <- glmer(Y ~ Treatment*Year + Skin + Age + Exposure + (1 | ID), data = skin, family = "poisson")
    skin2_sum = summary(skin2_glmm)
    skin2_sum$coefficients %>% knitr::kable()
